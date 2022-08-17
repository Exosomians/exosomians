options(stringsAsFactors = F)
options(scipen = 999)

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(stringr)
library(limma)
library(rtracklayer)

# metadata.tsv is downloaded from ENCODE website:
# https://www.encodeproject.org/metadata/?type=Experiment&status=released&assay_slims=RNA+binding&target.investigated_as=RNA+binding+protein
encode = read.delim('Data/Extreme/metadata.tsv', sep = '\t')
# View(encode)

# cleanings and getting human proteins only
encode = encode[grepl('human', encode$Experiment.target), ]
encode = encode[encode$Experiment.target!='Input library control-human',]
encode = encode[encode$Experiment.target!='No protein target control-human',]
encode = encode[!grepl('Control', encode$Experiment.target), ]
encode$Experiment.target = lapply(encode$Experiment.target, function(inp) {
  strsplit2(inp, '-| ')[1, 1]
})
encode = encode[grepl('signal', encode$Output.type), ]
# remove redundant columns
encode = encode[, -c(6:18)]
allProtList = unique(unlist(encode$Experiment.target))
encode = encode[encode$Assembly=='GRCh38', ]
hg38ProtList = unique(unlist(encode$Experiment.target))

# TODO: missed proteins
todoListProt = setdiff(allProtList, hg38ProtList)
# dim(encode)

# copy all files to coresponding directories
apply(encode, 1, function(encEntry) {
  targetProtein = encEntry['Experiment.target']
  spDir = paste('Data/Extreme/encode/bigWigAnalysis/files', targetProtein, sep = '/')
  fromDir = '/media/pgdrive/sharif/exosomians/data/ENCODE/files/bigwigs'
  if (!dir.exists(spDir)) {
    dir.create(spDir)
  }
  strandness = 'neg'
  if (grepl('plus', encEntry['Output.type'])) {
    strandness = 'pos'
  }
  fileAccession = encEntry['File.accession']
  fromDir = paste(fromDir, fileAccession, sep = '/')
  fromDir = paste(fromDir, 'bigWig', sep = '.')
  spDir = paste(spDir, strandness, sep = '/')
  if (!dir.exists(spDir)) {
    dir.create(spDir)
  }
  spDir = paste(spDir, fileAccession, sep = '/')
  spDir = paste(spDir, 'bigWig', sep = '.')
  file.symlink(from = fromDir, to = spDir)
})

preDir = 'Data/Extreme/encode/bigWigAnalysis'
baseDir = paste(preDir, 'files', sep = '/')
proteinDirectories = list.dirs(path = baseDir, full.names = FALSE, recursive = FALSE)
lapply(proteinDirectories, function(protDir) {
  specDir = paste(baseDir, protDir, 'pos', sep = '/')
  if (dir.exists(specDir)) {
    bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = TRUE, recursive = FALSE)
    lapply(bigwigFiles, function(bigwigInp){
      testbw = import.bw(bigwigInp)
      rm(testbw)
    })
  }
  specDir = paste(baseDir, protDir, 'neg', sep = '/')
  if (dir.exists(specDir)) {
    bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = TRUE, recursive = FALSE)
    lapply(bigwigFiles, function(bigwigInp){
      testbw = import.bw(bigwigInp)
      rm(testbw)
    })
  }
})

# import and create IC and EV extreme ranges
icExtPredictions = read.csv('Data/Extreme/ExoCNN.ic.extreme.90.probabilities.csv')
evExtPredictions = read.csv('Data/Extreme/ExoCNN.ev.extreme.90.probabilities.csv')

# head(icExtPredictions)
# View(icExtPredictions)
icExtPredictions = data.frame(str_split_fixed(icExtPredictions$id, '_', 5))
icExtPredictions = icExtPredictions[, c(1, 2, 3, 4)]
colnames(icExtPredictions) = c('chr', 'start', 'end', 'strand')
icExtPredictions$start = as.numeric(icExtPredictions$start)
icExtPredictions$end = as.numeric(icExtPredictions$end)

evExtPredictions = data.frame(str_split_fixed(evExtPredictions$id, '_', 5))
evExtPredictions = evExtPredictions[, c(1, 2, 3, 4)]
colnames(evExtPredictions) = c('chr', 'start', 'end', 'strand')
evExtPredictions$start = as.numeric(evExtPredictions$start)
evExtPredictions$end = as.numeric(evExtPredictions$end)

ICranges = GRanges(seqnames = icExtPredictions$chr, ranges = IRanges(start = icExtPredictions$start, end = icExtPredictions$end), strand = icExtPredictions$strand)
EVranges = GRanges(seqnames = evExtPredictions$chr, ranges = IRanges(start = evExtPredictions$start, end = evExtPredictions$end), strand = evExtPredictions$strand)

proteinDirectories = list.dirs(path = baseDir, full.names = FALSE, recursive = FALSE)
proteinStats = mclapply(proteinDirectories, function(protDir) {
  pos_merged_output = NULL
  specDir = paste(baseDir, protDir, 'pos', sep = '/')
  if (dir.exists(specDir)) {
    bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = FALSE, recursive = FALSE)
    if (length(bigwigFiles) > 1) {
      bigwigMergeTool = paste(preDir, 'bigWigMerge', sep = '/')
      bigwigMergedOutputName = paste(protDir, 'pos_unsorted_merged_output.bed', sep = '_')
      bigwidDirFiles = paste(specDir, '*.bigWig', sep = '/')
      bigwigMergeCommand = paste(bigwigMergeTool, bigwidDirFiles, bigwigMergedOutputName, sep = ' ')
      system(bigwigMergeCommand)
      sortTool = 'sort -k1,1 -k2,2n'
      sortedMergedOutputName = paste(protDir, 'pos_sorted_merged_output.bed', sep = '_')
      sortCommand = paste(sortTool, bigwigMergedOutputName, '>', sortedMergedOutputName, sep = ' ')
      system(sortCommand)
      bedToBigwigTool = paste(preDir, 'bedGraphToBigWig', sep = '/')
      # chrom size list file downloaded from:
      # https://www.encodeproject.org/files/GRCh38_EBV.chrom.sizes/
      chromSizes = paste(preDir, 'GRCh38_EBV.chrom.sizes', sep = '/')
      finalBigwigOutputName = paste(protDir, 'pos_merged_output.bigWig', sep = '_')
      bedToBigwigCommand = paste(bedToBigwigTool, sortedMergedOutputName, chromSizes, finalBigwigOutputName, sep = ' ')
      system(bedToBigwigCommand)
      pos_merged_output = import.bw(finalBigwigOutputName)
      # cleanings
      removeSortedOutputCommand = paste('rm', sortedMergedOutputName, sep = ' ')
      system(removeSortedOutputCommand)
      removeUnsortedBedOutputCommand = paste('rm', bigwigMergedOutputName, sep = ' ')
      system(removeUnsortedBedOutputCommand)
      removeFinalBigwigOutputCommand = paste('rm', finalBigwigOutputName, sep = ' ')
      system(removeFinalBigwigOutputCommand)
    } else {
      bigwidDirFiles = paste(specDir, bigwigFiles, sep = '/')
      pos_merged_output = import.bw(bigwidDirFiles)
    }
    pos_merged_output = as.data.frame(pos_merged_output)
    pos_merged_output$strand = NULL
    pos_merged_output$strand = '+'
  }
  neg_merged_output = NULL
  specDir = paste(baseDir, protDir, 'neg', sep = '/')
  if (dir.exists(specDir)) {
    bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = TRUE, recursive = FALSE)
    if (length(bigwigFiles) > 1) {
      lapply(bigwigFiles, function(bigwigInp) {
        testbw = import.bw(bigwigInp)
        testbw$score = abs(testbw$score)
        newName = paste(bigwigInp, 'adj.bw', sep = '_')
        export.bw(testbw, newName, format='bigWig')
        # cleanings
        removePrevBigWigCommand = paste('rm', bigwigInp, sep = ' ')
        system(removePrevBigWigCommand)
      })
      bigwigMergeTool = paste(preDir, 'bigWigMerge', sep = '/')
      bigwigMergedOutputName = paste(protDir, 'neg_unsorted_merged_output.bed', sep = '_')
      specDirFiles = paste(specDir, '*_adj.bw', sep = '/')
      bigwigMergeCommand = paste(bigwigMergeTool, specDirFiles, bigwigMergedOutputName, sep = ' ')
      system(bigwigMergeCommand)
      sortTool = 'sort -k1,1 -k2,2n'
      sortedMergedOutputName = paste(protDir, 'neg_sorted_merged_output.bed', sep = '_')
      sortCommand = paste(sortTool, bigwigMergedOutputName, '>', sortedMergedOutputName, sep = ' ')
      system(sortCommand)
      bedToBigwigTool = paste(preDir, 'bedGraphToBigWig', sep = '/')
      chromSizes = paste(preDir, 'GRCh38_EBV.chrom.sizes', sep = '/')
      finalBigwigOutputName = paste(protDir, 'neg_merged_output.bigWig', sep = '_')
      bedToBigwigCommand = paste(bedToBigwigTool, sortedMergedOutputName, chromSizes, finalBigwigOutputName, sep = ' ')
      system(bedToBigwigCommand)
      neg_merged_output = import.bw(finalBigwigOutputName)
      # cleanings
      removeSortedOutputCommand = paste('rm', sortedMergedOutputName, sep = ' ')
      system(removeSortedOutputCommand)
      removeUnsortedBedOutputCommand = paste('rm', bigwigMergedOutputName, sep = ' ')
      system(removeUnsortedBedOutputCommand)
      removeFinalBigwigOutputCommand = paste('rm', finalBigwigOutputName, sep = ' ')
      system(removeFinalBigwigOutputCommand)
      removeSpecDirFilesCommand = paste('rm', specDirFiles, sep = ' ')
      system(removeSpecDirFilesCommand)
    } else {
      neg_merged_output = import.bw(bigwigFiles)
      neg_merged_output$score = abs(neg_merged_output$score)
    }
    neg_merged_output = as.data.frame(neg_merged_output)
    neg_merged_output$strand = NULL
    neg_merged_output$strand = '-'
  }
  merged_output = rbind(neg_merged_output, pos_merged_output)
  merged_ranges = GRanges(seqnames = merged_output$seqnames, ranges = IRanges(start = merged_output$start, end = merged_output$end), strand = merged_output$strand)
  
  # findOverlaps
  merged_output$inEV = overlapsAny(merged_ranges, EVranges)
  merged_output$inIC = overlapsAny(merged_ranges, ICranges)
  
  icRegs = merged_output[merged_output$inIC, ]
  icMean = weighted.mean(x = icRegs$score, w = icRegs$width)
  evRegs = merged_output[merged_output$inEV, ]
  evMean = weighted.mean(x = evRegs$score, w = evRegs$width)
  totMean = weighted.mean(x = merged_output$score, w = merged_output$width)
  threshold = quantile(x = merged_output$score, probs = 0.3)
  diffr = evMean - icMean
  result = paste(evMean, icMean, totMean, diffr, threshold, sep = '_')
  if (diffr > threshold && diffr > totMean) {
    result = paste(result, 'EV', sep = '_')
    return(result)
  } else if (abs(diffr) > threshold && abs(diffr) > totMean) {
    result = paste(result, 'IC', sep = '_')
    return(result)
  } else {
    result = paste(result, 'None', sep = '_')
    return(result)
  }
}, mc.cores = 7)

finalResult = data.frame(Protein=proteinDirectories,
                         Stat=unlist(proteinStats))
tmpFinalResult = data.frame(str_split_fixed(finalResult$Stat, '_', 6))
finalResult = cbind(finalResult$Protein, tmpFinalResult)
colnames(finalResult) = c('Protein', 'evMean', 'icMean', 'totMean', 'difference', 'threshold', 'verdict')
saveRDS(finalResult, 'finalResult.rds')

