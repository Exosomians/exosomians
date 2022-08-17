setwd('/media/pgdrive/sharif/exosomians/data/bams')

options(stringsAsFactors = F)
# source('Codes/Functions.R')
# Initialize()
library(parallel)

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 <- BSgenome.Hsapiens.UCSC.hg38

#### Extracting features using annotation files ####
IC_BED = 'ic.final.bed'
EV_BED = 'ev.final.bed'

bedFilePath = EV_BED
label = 'ev'
genome = hg38


evDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = EV_BED,
                                                          label = 'ev',
                                                          genome = hg38)

icDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = IC_BED,
                                                          label = 'ic',
                                                          genome = hg38)

# dim(evDesignMat)
# dim(icDesignMat)


designMat1 = Preprocess.Downsample(icDesignMat,
                                   evDesignMat,
                                   seed = 1)
designMat2 = Preprocess.Downsample(icDesignMat,
                                   evDesignMat,
                                   seed = 2)
designMat3 = Preprocess.Downsample(icDesignMat,
                                   evDesignMat,
                                   seed = 3)


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
  smRNAsSeq <- getSeq(hg38, smRNAsRange)
  
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
