setwd('/media/pgdrive/sharif/exosomians/data/bams')

options(stringsAsFactors = F)
# source('Codes/Functions.R')
# Initialize()
library(parallel)

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

#### Extracting features using annotation files ####
EV_LABEL = 'ev'
EV_BED = 'ev.final.bed'
EV_FASTA = 'ev.final.fasta'


evDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = EV_BED,
                                                          label = EV_LABEL,
                                                          genome = hg38)

Preprocess.Write2Fasta(designMat = evDesignMat,
                       outputFile = EV_FASTA)

system("bash -xv main-vienna.sh 'ev.final.fasta'")
evDesignMat = Preprocess.AddDotBracketSeq(evDesignMat,
                                          dotBracketFasta = gsub('.fasta',
                                                                 '.dotbracket.fasta',
                                                                 EV_FASTA))

write.csv(evDesignMat, file = paste0(EV_LABEL, '_DesignMatrix_SeqPlusDotBracket.csv'),
          quote = F, row.names = F)


Preprocess.SequenceExtractionAndRefinements = function(bedFilePath, label, genome)
{
  bedFile <- read.delim(bedFilePath, header = F)
  colnames(bedFile) <- c('seqnames', 'start', 'end', 'name', 'score', 'strand')
  bedFile$name <- paste0(bedFile$seqnames, '_', bedFile$start, '_', bedFile$end, '_', bedFile$strand, '_', label)
  bedFile <- subset(bedFile, seqnames %in% paste0('chr', c(1:22,'X', 'Y')))
  
  smRNAsRange = GRangesForBSGenome(genome = 'hg38',
                                   chrom = bedFile$seqnames,
                                   ranges = IRanges(start = bedFile$start,
                                                    end = bedFile$end),
                                   strand = bedFile$strand)
  
  #### Sequence Extraction
  smRNAsSeq = getSeq(hg38, smRNAsRange)
  smRNAsSeq = RNAStringSet(smRNAsSeq)
  
  #### Sequence Refinements
  
  ### Correcting sequences based on their strand
  positiveStrandsSeqsIndex = which(bedFile$strand=='+')
  negativeStrandsSeqsIndex = which(bedFile$strand=='-')
  smRNAsSeq[positiveStrandsSeqsIndex] = complement(smRNAsSeq[positiveStrandsSeqsIndex])
  smRNAsSeq[negativeStrandsSeqsIndex] = reverseComplement(smRNAsSeq[negativeStrandsSeqsIndex])
  
  designMat = data.frame(id = bedFile$name,
                         chr = bedFile$seqnames,
                         strand = bedFile$strand,
                         seq = as.character(smRNAsSeq),
                         length = width(smRNAsSeq),
                         label = label)
  
  ##### Removing sequences that have N.
  sequencesContainsN = mclapply(strsplit(designMat$seq, ''),
                                function(eachSequence) 'N'%in%eachSequence,
                                mc.cores = detectCores()-2)
  sequencesContainsN = which(unlist(sequencesContainsN))
  
  designMat = designMat[-sequencesContainsN, ]
  
  
  if(any(duplicated(designMat$id)))
    print('CAUTION! Duplicated Ids Exist!')
  
  designMat
}


Preprocess.Downsample = function(icDesignMat, evDesignMat, seed)
{
  set.seed(seed)
  
  numOfSamples = min(nrow(icDesignMat), nrow(evDesignMat))
  
  icSamplesIndexToElicit = sample(seq(nrow(icDesignMat)), numOfSamples)
  evSamplesIndexToElicit = sample(seq(nrow(evDesignMat)), numOfSamples)
  
  rbind(icDesignMat[icSamplesIndexToElicit, ],
        evDesignMat[evSamplesIndexToElicit, ])
}


Preprocess.Write2Fasta = function(designMat, outputFile)
{
  seqinr::write.fasta(sequences = strsplit(designMat$seq, ''),
                      names = designMat$id,
                      file.out = outputFile)
}

Preprocess.AddDotBracketSeq = function(designMat, dotBracketFasta)
{
  secStr = readBStringSet(dotBracketFasta)
  secStr = substr(secStr, 1, designMat$length)
  designMat$dotbracket = secStr
  designMat
}
