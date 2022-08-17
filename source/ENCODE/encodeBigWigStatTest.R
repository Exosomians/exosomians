options(stringsAsFactors = F)
# for scientific representations:
# options(scipen = 50)
options(scipen = 999)

library(BSgenome.Hsapiens.UCSC.hg38)
hg38 = BSgenome.Hsapiens.UCSC.hg38
library(stringr)
library(limma)
library(rtracklayer)
library(data.table)

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

preDir = 'Data/Extreme/encode/bigWigAnalysis'
baseDir = paste(preDir, 'files', sep = '/')
outDir = 'Data/Extreme/encode/bigWigOutputs'

finalResult = readRDS('finalResult.rds')
finalResult$difference = as.numeric(finalResult$difference)
finalResult = finalResult[order(finalResult$difference, decreasing = T), ]

pVals = mclapply(finalResult$Protein, function(protDir) {
  mergedOutputName = paste(protDir, 'peaked', 'rds', sep = '.')
  mergedOutputPath = paste(outDir, mergedOutputName, sep = '/')
  if (!file.exists(mergedOutputPath)) {
    return(-1)
  }
  merged_output = readRDS(mergedOutputPath)
  # View(merged_output)
  # dim(merged_output)
  # summary(merged_output$score)
  # quantile(merged_output$score, seq(0.9, 1, 0.01))
  # merged_output = merged_output[merged_output$score>=median(merged_output$score), ]
  merged_output = merged_output[merged_output$isPeak, ]
  # icRegs = merged_output[merged_output$inIC, ]
  # nonicRegs = merged_output[!merged_output$inIC, ]
  # evRegs = nonicRegs[nonicRegs$inEV, ]
  # merged_output = rbind(icRegs, evRegs)
  icRegs = merged_output[merged_output$inIC, ]
  nonicRegs = merged_output[!merged_output$inIC, ]
  nonicRegs$score = 0
  evRegs = merged_output[merged_output$inEV, ]
  nonevRegs = merged_output[!merged_output$inEV, ]
  nonevRegs$score = 0
  icScores = rep(icRegs$score, icRegs$width)
  icScores = c(icScores, rep(nonicRegs$score, nonicRegs$width))
  evScores = rep(evRegs$score, evRegs$width)
  evScores = c(evScores, rep(nonevRegs$score, nonevRegs$width))
  # quantile(icScores, seq(0.9995, 1, 0.00005))
  # quantile(evScores, seq(0.9995, 1, 0.00005))
  wilcoxTest = wilcox.test(x=evScores, y=icScores, alternative = 'greater', paired = F)
  wilcoxTest$p.value
}, mc.cores = 10)

pVals = lapply(pVals, function(x) ifelse(is.null(x), NA, x))
pVals = unlist(pVals)
# table(pVals)
finalResult$pVals = pVals
finalResult = finalResult[finalResult$pVals>=0,]
finalResult$adjPvals = p.adjust(finalResult$pVals, method = 'BH')
finalResult = finalResult[finalResult$adjPvals<=0.05,]
finalResult = finalResult[order(finalResult$adjPvals), ]
finalResult = finalResult[, c(1, 8, 9)]
saveRDS(finalResult, 'finalResultPeakedAdjPval.rds')

# saveRDS(pVals, 'medianPvalues.rds')
# pVals = readRDS('meanPvaluesSec.rds')

# metadata.tsv is downloaded from ENCODE website:
# https://www.encodeproject.org/metadata/?type=Experiment&status=released&assay_slims=RNA+binding&target.investigated_as=RNA+binding+protein
encode = read.delim('Data/Extreme/metadata.tsv', sep = '\t')
# View(encode)

# extracting human proteins only
encode = encode[grepl('human', encode$Experiment.target), ]
encode = encode[encode$Experiment.target!='Input library control-human',]
encode = encode[encode$Experiment.target!='No protein target control-human',]
encode = encode[!grepl('Control', encode$Experiment.target), ]
encode$Experiment.target = lapply(encode$Experiment.target, function(inp) {
  strsplit2(inp, '-| ')[1, 1]
})

encodePeaks = encode[grepl('peaks', encode$Output.type), ]
encodePeaks = encodePeaks[!grepl('bigBed', encodePeaks$File.format), ]
encodePeaks = encodePeaks[encodePeaks$Assembly=='GRCh38', ]
encodePeaks = encodePeaks[encodePeaks$Biological.replicate.s.=='1, 2', ]
# dim(encodePeaks)
# View(encodePeaks)

# copy all files to coresponding directories
bedsBaseDir = paste(baseDir, 'beds', sep = '/')
apply(encodePeaks, 1, function(encEntry) {
  targetProtein = encEntry['Experiment.target']
  spDir = paste(bedsBaseDir, targetProtein, sep = '/')
  fromDir = '/media/pgdrive/sharif/exosomians/data/ENCODE/files/beds'
  if (!dir.exists(spDir)) {
    dir.create(spDir)
  }
  fileAccession = encEntry['File.accession']
  fromDir = paste(fromDir, fileAccession, sep = '/')
  fromDir = paste(fromDir, 'bed.gz', sep = '.')
  spDir = paste(spDir, fileAccession, sep = '/')
  spDir = paste(spDir, 'bed.gz', sep = '.')
  file.copy(from = fromDir, to = spDir)
})

# to unzip and check whether there exist any corrupted files
proteinDirectories = list.dirs(path = bedsBaseDir, full.names = FALSE, recursive = FALSE)
lapply(proteinDirectories, function(protDir) {
  specDir = paste(bedsBaseDir, protDir, sep = '/')
  if (dir.exists(specDir)) {
    bedFiles = list.files(path = specDir, pattern = "\\.bed.gz$", full.names = TRUE, recursive = FALSE)
    lapply(bedFiles, function(bedInp) {
      unzipCommand = paste('gunzip', bedInp, sep = ' ')
      system(unzipCommand)
    })
  }
})

# to merge bed files of a specific protein
proteinDirectories = list.dirs(path = bedsBaseDir, full.names = FALSE, recursive = FALSE)
lapply(proteinDirectories, function(protDir) {
  specDir = paste(bedsBaseDir, protDir, sep = '/')
  if (dir.exists(specDir)) {
    bedFiles = list.files(path = specDir, pattern = "\\.bed$", full.names = TRUE, recursive = FALSE)
    if (length(bedFiles) > 1) {
      # merging
      allBeds = paste(specDir, '*.bed', sep = '/')
      allBedsFilename = paste('all', protDir, 'bed', sep = '.')
      allBedsFilename = paste(specDir, allBedsFilename, sep = '/')
      mergeCommand = paste('cat', allBeds, '>', allBedsFilename, sep = ' ')
      system(mergeCommand)
      lapply(bedFiles, function(bedInp) {
        rmCommand = paste('rm', bedInp, sep = ' ')
        system(rmCommand)
      })
    }
  }
})

# these proteins do not have any peaks file
uncoveredProts = finalResult$Protein[!finalResult$Protein%in%proteinDirectories]

# adding isPeak column to previously merged regions
mclapply(proteinDirectories, function(protDir) {
  mergedOutputName = paste(protDir, 'rds', sep = '.')
  mergedOutputPath = paste(outDir, mergedOutputName, sep = '/')
  merged_output = readRDS(mergedOutputPath)
  merged_ranges = GRanges(seqnames = merged_output$seqnames, ranges = IRanges(start = merged_output$start, end = merged_output$end), strand = merged_output$strand)
  
  specDir = paste(bedsBaseDir, protDir, sep = '/')
  bedFile = list.files(path = specDir, pattern = "\\.bed$", full.names = TRUE, recursive = FALSE)
  peakRegions = fread(bedFile)
  peakRegions = peakRegions[, c(1, 2, 3, 6)]
  colnames(peakRegions) = c('chr', 'start', 'end', 'strand')
  peakRegions$start = as.numeric(peakRegions$start)
  peakRegions$end = as.numeric(peakRegions$end)
  peakRanges = GRanges(seqnames = peakRegions$chr, ranges = IRanges(start = peakRegions$start, end = peakRegions$end), strand = peakRegions$strand)
  
  # findOverlaps
  merged_output$isPeak = overlapsAny(merged_ranges, peakRanges)
  newMergedOutputName = paste(protDir, 'peaked', 'rds', sep = '.')
  newMergedOutputPath = paste(outDir, newMergedOutputName, sep = '/')
  saveRDS(merged_output, newMergedOutputPath)
  rmCommand = paste('rm', mergedOutputPath)
  system(rmCommand)
}, mc.cores = 10)

# extracting extreme EV sequences that have overlap with peaked regions
finalResult = readRDS('finalResultPeakedAdjPval.rds')
finalResult = finalResult[finalResult$adjPvals<=0.05,]
mclapply(finalResult$Protein, function(protDir) {
  # protDir = 'NSUN2'
  mergedOutputName = paste(protDir, 'peaked', 'rds', sep = '.')
  mergedOutputPath = paste(outDir, mergedOutputName, sep = '/')
  merged_output = readRDS(mergedOutputPath)
  merged_output = merged_output[merged_output$inEV, ]
  merged_output = merged_output[merged_output$isPeak, ]
  merged_output = merged_output[, seq(1, 6)]
  outputBedName = paste(protDir, 'ev', 'peaked', 'bed', sep = '.')
  outputBedPath = paste(outDir, outputBedName, sep = '/')
  fwrite(merged_output, outputBedPath, col.names = F, sep = '\t')
  outputBedSortedName = paste(protDir, 'ev', 'peaked', 'sorted', 'bed', sep = '.')
  outputBedSortedPath = paste(outDir, outputBedSortedName, sep = '/')
  sortTool = 'sort -k1,1 -k2,2n'
  sortCommand = paste(sortTool, outputBedPath, '>', outputBedSortedPath, sep = ' ')
  system(sortCommand)
  posOutputBedName = paste(protDir, 'ev', 'peaked', 'pos', 'bed', sep = '.')
  posOutputBedPath = paste(outDir, posOutputBedName, sep = '/')
  negOutputBedName = paste(protDir, 'ev', 'peaked', 'neg', 'bed', sep = '.')
  negOutputBedPath = paste(outDir, negOutputBedName, sep = '/')
  mergePosBedCommand = paste('bedtools merge -i', outputBedSortedPath, '-d 1 -S +', '>', posOutputBedPath, sep = ' ')
  system(mergePosBedCommand)
  posOutput = fread(posOutputBedPath)
  if(nrow(posOutput) > 0) {
    colnames(posOutput) = c('chr', 'start', 'end', 'strand')
    posOutput$name = '.'
    posOutput$score = 0
    posOutput = posOutput[, c(1, 2, 3, 5, 6, 4)]
  } else {
    posOutput = NULL
  }
  mergeNegBedCommand = paste('bedtools merge -i', outputBedSortedPath, '-d 1 -S -', '>', negOutputBedPath, sep = ' ')
  system(mergeNegBedCommand)
  negOutput = fread(negOutputBedPath)
  if(nrow(negOutput) > 0) {
    colnames(negOutput) = c('chr', 'start', 'end', 'strand')
    negOutput$name = '.'
    negOutput$score = 0
    negOutput = negOutput[, c(1, 2, 3, 5, 6, 4)]
  } else {
    negOutput = NULL
  }
  allMergedOutput = rbind(negOutput, posOutput)
  allMergedOutputRanges = GRanges(seqnames = allMergedOutput$chr, ranges = IRanges(start = allMergedOutput$start, end = allMergedOutput$end), strand = allMergedOutput$strand)
  
  # get a copy of evExtPredictions
  evpredicns = evExtPredictions
  
  # findOverlaps
  evpredicns$isPeak = overlapsAny(EVranges, allMergedOutputRanges)
  evpredicns = evpredicns[evpredicns$isPeak, ]
  # extracting extreme EV sequence IDs
  if (nrow(evpredicns) > 0) {
    evpredicns$id = paste(evpredicns$chr, evpredicns$start, evpredicns$end, evpredicns$strand, sep = '_')
    evpredicns$Protein = protDir
    evpredicns = evpredicns[, c('id', 'Protein')]
  } else {
    evpredicns = NULL
  }
  outputRegName = paste(protDir, 'ev', 'peaked', 'regions', sep = '.')
  outputRegPath = paste(outDir, outputRegName, sep = '/')
  fwrite(evpredicns, outputRegPath, sep = '\t')
  rmCommand = paste('rm', outputBedPath)
  system(rmCommand)
  rmCommand = paste('rm', posOutputBedPath)
  system(rmCommand)
  rmCommand = paste('rm', negOutputBedPath)
  system(rmCommand)
  rmCommand = paste('rm', outputBedSortedPath)
  system(rmCommand)
}, mc.cores = 5)

allEVPeakedRegions = NULL
for(prot in finalResult$Protein) {
  # print(prot)
  outputRegName = paste(prot, 'ev', 'peaked', 'regions', sep = '.')
  outputRegPath = paste(outDir, outputRegName, sep = '/')
  protRegs = fread(outputRegPath)
  allEVPeakedRegions = rbind(allEVPeakedRegions, protRegs)
}
outputRegPath = paste(outDir, 'allEvProtsOverlapEvExtreme', sep = '/')
# fwrite(allEVPeakedRegions, outputRegPath)
allEVPeakedRegions = fread(outputRegPath)
length(unique(allEVPeakedRegions$Protein))

# extracting extreme EV and IC sequences that have overlap with peaked regions,
# and their average coverage per each
finalResult = readRDS('finalResultPeakedAdjPval.rds')
finalResult = finalResult[finalResult$adjPvals<=0.05,]
mclapply(finalResult$Protein, function(protDir) {
  # protDir = 'NSUN2'
  mergedOutputName = paste(protDir, 'peaked', 'rds', sep = '.')
  mergedOutputPath = paste(outDir, mergedOutputName, sep = '/')
  merged_output = readRDS(mergedOutputPath)
  merged_output = merged_output[merged_output$isPeak, ]
  merged_output_ev = merged_output[merged_output$inEV, ]
  merged_output_ic = merged_output[merged_output$inIC, ]
  
  icMean = weighted.mean(x = merged_output_ic$score, w = merged_output_ic$width)
  if (is.nan(icMean)) {
    icMean = 0
  }
  evMean = weighted.mean(x = merged_output_ev$score, w = merged_output_ev$width)
  if (is.nan(evMean)) {
    evMean = 0
  }
  
  mergedOutputRangesEV = GRanges(seqnames = merged_output_ev$seqnames, ranges = IRanges(start = merged_output_ev$start, end = merged_output_ev$end), strand = merged_output_ev$strand)
  mergedOutputRangesIC = GRanges(seqnames = merged_output_ic$seqnames, ranges = IRanges(start = merged_output_ic$start, end = merged_output_ic$end), strand = merged_output_ic$strand)
  # get a copy of evExtPredictions and icExtPredictions
  evpredicns = evExtPredictions
  icpredicns = icExtPredictions
  # findOverlaps
  evpredicns$isPeak = overlapsAny(EVranges, mergedOutputRangesEV)
  icpredicns$isPeak = overlapsAny(ICranges, mergedOutputRangesIC)
  evpredicns = evpredicns[evpredicns$isPeak, ]
  icpredicns = icpredicns[icpredicns$isPeak, ]
  # extracting extreme EV sequence IDs
  if (nrow(evpredicns) > 0) {
    evpredicns$id = paste(evpredicns$chr, evpredicns$start, evpredicns$end, evpredicns$strand, sep = '_')
    evpredicns$seq = as.character(RNAStringSet(getSeq(hg38, names = evpredicns$chr, start = evpredicns$start, end = evpredicns$end, strand = evpredicns$strand)))
    evpredicns$type = 'EV'
    evpredicns$Protein = protDir
    evpredicns$evMean = evMean
    evpredicns$icMean = icMean
    evpredicns = evpredicns[, c('id', 'seq', 'type', 'Protein', 'evMean', 'icMean')]
  } else {
    evpredicns = NULL
  }
  # extracting extreme IC sequence IDs
  if (nrow(icpredicns) > 0) {
    icpredicns$id = paste(icpredicns$chr, icpredicns$start, icpredicns$end, icpredicns$strand, sep = '_')
    icpredicns$seq = as.character(RNAStringSet(getSeq(hg38, names = icpredicns$chr, start = icpredicns$start, end = icpredicns$end, strand = icpredicns$strand)))
    icpredicns$type = 'IC'
    icpredicns$Protein = protDir
    icpredicns$evMean = evMean
    icpredicns$icMean = icMean
    icpredicns = icpredicns[, c('id', 'seq', 'type', 'Protein', 'evMean', 'icMean')]
  } else {
    icpredicns = NULL
  }
  EVandICmerged = rbind(evpredicns, icpredicns)
  outputRegName = paste(protDir, 'EVandIC', 'peaked', 'regions', sep = '.')
  outputRegPath = paste(outDir, outputRegName, sep = '/')
  fwrite(EVandICmerged, outputRegPath, sep = '\t')
}, mc.cores = 5)

allEVandICPeakedRegions = NULL
for(prot in finalResult$Protein) {
  # print(prot)
  outputRegName = paste(prot, 'EVandIC', 'peaked', 'regions', sep = '.')
  outputRegPath = paste(outDir, outputRegName, sep = '/')
  protRegs = fread(outputRegPath)
  allEVandICPeakedRegions = rbind(allEVandICPeakedRegions, protRegs)
}
# View(allEVandICPeakedRegions)

outputRegPath = paste(outDir, 'EvProtsOverlapEvandIcExtreme', sep = '/')
# fwrite(allEVandICPeakedRegions, outputRegPath)
allEVandICPeakedRegions = as.data.table(allEVandICPeakedRegions)
class(allEVandICPeakedRegions)
dim(allEVandICPeakedRegions)
evTypedPeakedRegions = allEVandICPeakedRegions[allEVandICPeakedRegions$type=='EV',]
icTypedPeakedRegions = allEVandICPeakedRegions[allEVandICPeakedRegions$type=='IC',]
evTypedPeakedRegions = evTypedPeakedRegions[, list(EVseqsNo=length(unique(id)), EVseqIDs=toString(unique(id)), EVseqs=toString(unique(seq)), EVmeanPeaks=max(evMean), ICmeanPeaks=max(icMean)), by=list(Protein)]
icTypedPeakedRegions = icTypedPeakedRegions[, list(ICseqsNo=length(unique(id)), ICseqIDs=toString(unique(id)), ICseqs=toString(unique(seq))), by=list(Protein)]

allTypePeakedRegions = merge(x=evTypedPeakedRegions, y=icTypedPeakedRegions, by='Protein', all=T)
finalResult = finalResult[, c('Protein', 'pVals', 'adjPvals')]
allTypePeakedRegions = merge(x=allTypePeakedRegions, y=finalResult, by='Protein', all=T)
allTypePeakedRegions$adjPvals = as.numeric(allTypePeakedRegions$adjPvals)
allTypePeakedRegions = allTypePeakedRegions[order(allTypePeakedRegions$adjPvals), ]
# View(allTypePeakedRegions)
# saveRDS(allTypePeakedRegions, 'allTypePeakedRegions')
write.table(allTypePeakedRegions, file='allTypePeakedRegions.tsv', row.names=FALSE, sep='\t', quote=F)

allTypePeakedRegionsPreview = allTypePeakedRegions[, c('Protein', 'EVseqsNo', 'ICseqsNo', 'EVmeanPeaks', 'ICmeanPeaks', 'adjPvals')]
allTypePeakedRegionsPreview$adjPvals = as.numeric(allTypePeakedRegionsPreview$adjPvals)
allTypePeakedRegionsPreview = allTypePeakedRegionsPreview[order(allTypePeakedRegionsPreview$adjPvals), ]
# View(allTypePeakedRegionsPreview)
# saveRDS(allTypePeakedRegionsPreview, 'allTypePeakedRegionsPreview')





