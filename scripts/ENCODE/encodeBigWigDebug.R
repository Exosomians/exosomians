protDir = 'IGF2BP3'

print(protDir)
pos_merged_output = NULL
specDir = paste(baseDir, protDir, 'pos', sep = '/')
if (dir.exists(specDir)) {
  # setwd(specDir)
  bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = FALSE, recursive = FALSE)
  if (length(bigwigFiles) > 1) {
    bigwigMergeTool = paste(preDir, 'bigWigMerge', sep = '/')
    bigwigMergedOutputName = paste(protDir, 'pos_unsorted_merged_output.bed', sep = '_')
    bigwidDirFiles = paste(specDir, '*.bigWig', sep = '/')
    bigwigMergeCommand = paste(bigwigMergeTool, bigwidDirFiles, bigwigMergedOutputName, sep = ' ')
    print(bigwigMergeCommand)
    system(bigwigMergeCommand)
    sortTool = 'sort -k1,1 -k2,2n'
    sortedMergedOutputName = paste(protDir, 'pos_sorted_merged_output.bed', sep = '_')
    sortCommand = paste(sortTool, bigwigMergedOutputName, '>', sortedMergedOutputName, sep = ' ')
    print(sortCommand)
    system(sortCommand)
    bedToBigwigTool = paste(preDir, 'bedGraphToBigWig', sep = '/')
    chromSizes = paste(preDir, 'GRCh38_EBV.chrom.sizes', sep = '/')
    finalBigwigOutputName = paste(protDir, 'pos_merged_output.bigWig', sep = '_')
    bedToBigwigCommand = paste(bedToBigwigTool, sortedMergedOutputName, chromSizes, finalBigwigOutputName, sep = ' ')
    print(bedToBigwigCommand)
    system(bedToBigwigCommand)
    print(finalBigwigOutputName)
    pos_merged_output = import.bw(finalBigwigOutputName)
    # cleanings
    removeSortedOutputCommand = paste('rm', sortedMergedOutputName, sep = ' ')
    print(removeSortedOutputCommand)
    system(removeSortedOutputCommand)
    removeUnsortedBedOutputCommand = paste('rm', bigwigMergedOutputName, sep = ' ')
    print(removeUnsortedBedOutputCommand)
    system(removeUnsortedBedOutputCommand)
    removeFinalBigwigOutputCommand = paste('rm', finalBigwigOutputName, sep = ' ')
    print(removeFinalBigwigOutputCommand)
    system(removeFinalBigwigOutputCommand)
  } else {
    bigwidDirFiles = paste(specDir, bigwigFiles, sep = '/')
    print(bigwidDirFiles)
    pos_merged_output = import.bw(bigwidDirFiles)
  }
  pos_merged_output = as.data.frame(pos_merged_output)
  pos_merged_output$strand = NULL
  pos_merged_output$strand = '+'
}
# setwd(homeDir)
neg_merged_output = NULL
specDir = paste(baseDir, protDir, 'neg', sep = '/')
if (dir.exists(specDir)) {
  # setwd(specDir)
  bigwigFiles = list.files(path = specDir, pattern = "\\.bigWig$", full.names = TRUE, recursive = FALSE)
  if (length(bigwigFiles) > 1) {
    lapply(bigwigFiles, function(bigwigInp){
      print(bigwigInp)
      testbw = import.bw(bigwigInp)
      testbw$score = abs(testbw$score)
      newName = paste(bigwigInp, 'adj.bw', sep = '_')
      print('newName for new BigWig:')
      print(newName)
      export.bw(testbw, newName, format='bigWig')
      # cleanings
      removePrevBigWigCommand = paste('rm', bigwigInp, sep = ' ')
      print(removePrevBigWigCommand)
      system(removePrevBigWigCommand)
    })
    bigwigMergeTool = paste(preDir, 'bigWigMerge', sep = '/')
    bigwigMergedOutputName = paste(protDir, 'neg_unsorted_merged_output.bed', sep = '_')
    specDirFiles = paste(specDir, '*_adj.bw', sep = '/')
    bigwigMergeCommand = paste(bigwigMergeTool, specDirFiles, bigwigMergedOutputName, sep = ' ')
    print(bigwigMergeCommand)
    system(bigwigMergeCommand)
    sortTool = 'sort -k1,1 -k2,2n'
    sortedMergedOutputName = paste(protDir, 'neg_sorted_merged_output.bed', sep = '_')
    sortCommand = paste(sortTool, bigwigMergedOutputName, '>', sortedMergedOutputName, sep = ' ')
    print(sortCommand)
    system(sortCommand)
    bedToBigwigTool = paste(preDir, 'bedGraphToBigWig', sep = '/')
    chromSizes = paste(preDir, 'GRCh38_EBV.chrom.sizes', sep = '/')
    finalBigwigOutputName = paste(protDir, 'neg_merged_output.bigWig', sep = '_')
    bedToBigwigCommand = paste(bedToBigwigTool, sortedMergedOutputName, chromSizes, finalBigwigOutputName, sep = ' ')
    print(bedToBigwigCommand)
    system(bedToBigwigCommand)
    print(finalBigwigOutputName)
    neg_merged_output = import.bw(finalBigwigOutputName)
    # cleanings
    removeSortedOutputCommand = paste('rm', sortedMergedOutputName, sep = ' ')
    print(removeSortedOutputCommand)
    system(removeSortedOutputCommand)
    removeUnsortedBedOutputCommand = paste('rm', bigwigMergedOutputName, sep = ' ')
    print(removeUnsortedBedOutputCommand)
    system(removeUnsortedBedOutputCommand)
    removeFinalBigwigOutputCommand = paste('rm', finalBigwigOutputName, sep = ' ')
    print(removeFinalBigwigOutputCommand)
    system(removeFinalBigwigOutputCommand)
    removeSpecDirFilesCommand = paste('rm', specDirFiles, sep = ' ')
    print(removeSpecDirFilesCommand)
    system(removeSpecDirFilesCommand)
  } else {
    print(bigwigFiles)
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
  # return(result)
} else if (abs(diffr) > threshold && abs(diffr) > totMean) {
  result = paste(result, 'IC', sep = '_')
  # return(result)
} else {
  result = paste(result, 'None', sep = '_')
  # return(result)
}

head(merged_output)
table(merged_output$inEV)
table(merged_output$inIC)
head(result)
# View(merged_output[merged_output$inIC, ])
# View(merged_output[merged_output$inEV, ])
summary(merged_output[merged_output$inIC, ]$score)
summary(merged_output[merged_output$inEV, ]$score)

### Manipulating results
# View(head(merged_output))
# View(merged_output[merged_output$inEV==T, ])
inev = merged_output[merged_output$inEV==T, ]
nonev = merged_output[!merged_output$inEV, ]
nonev$score = 0
inic = merged_output[merged_output$inIC==T, ]
nonic = merged_output[!merged_output$inIC, ]
nonic$score = 0

icScores = rep(inic$score, inic$width)
icScores = c(icScores, rep(nonic$score, nonic$width))
evScores = rep(inev$score, inev$width)
evScores = c(evScores, rep(nonev$score, nonev$width))
sum(inic$width)
nsun2WilcoxTest = wilcox.test(x=evScores, y=icScores, alternative = 'greater')
IGF2BP3wilcoxTest = wilcox.test(x=evScores, y=icScores, alternative = 'greater')
IGF2BP3wilcoxTestPvalue = IGF2BP3wilcoxTest$p.value

weighted.mean(inev$score, inev$width)
sum(inev$score)
sum(inic$score)
sum(inev$width)
sum(inic$width)
View(inev)
summary(inev$score)
summary(inic$score)

inevAboveMedian = inev[inev$score > median(inev$score), ]
inicAboveMedian = inic[inic$score > median(inic$score), ]

weighted.mean(inevAboveMedian$score, inevAboveMedian$width)
weighted.mean(inicAboveMedian$score, inicAboveMedian$width)

quantile(inev$score, seq(0.75, 1, 0.05))
quantile(inic$score, seq(0.75, 1, 0.05))

# View(as.data.frame(EVranges))
# View(as.data.frame(merged_ranges))
# 
# View(pos_merged_output)
View(inev)
quantile(merged_output$score, seq(0.9995, 1, 0.00005))

