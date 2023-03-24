options(stringsAsFactors = F)
options(scipen = 999)

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(stringr)
library(limma)
library(rtracklayer)

preDir = 'Data/Extreme/encode/bigWigAnalysis'
baseDir = paste(preDir, 'files', sep = '/')
outDir = 'Data/Extreme/encode/bigWigOutputs'

finalResult = readRDS('finalResult.rds')
prots = finalResult$Protein
proteinStats = mclapply(prots, function(protDir) {
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
      lapply(bigwigFiles, function(bigwigInp){
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
  outname = paste(protDir, 'rds', sep = '.')
  outname = paste(outDir, outname, sep = '/')
  saveRDS(merged_output, outname)
}, mc.cores = 10)

