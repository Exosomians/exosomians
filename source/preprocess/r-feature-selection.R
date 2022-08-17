# setwd('/media/pgdrive/sharif/exosomians/data/bams')

options(stringsAsFactors = F)
source('r-functions.R')

library(parallel)
library(Biostrings)
library(limma)

finalResults = read.csv('final.result.mat.csv')
finalResults$X = NULL
# head(finalResults)

theLabels = strsplit2(finalResults$id, '_')[, 5]
# probEvYes = finalResults[theLabels=='ev', 'yes']
# probIcNo = finalResults[theLabels=='ic', 'no']

# hist(finalResults[theLabels=='ev', 'yes'])
# hist(finalResults[theLabels=='ic', 'no'])

# summary(finalResults[theLabels=='ev', 'yes'])
# summary(finalResults[theLabels=='ic', 'no'])

extremeEvsThreshold = quantile(finalResults[theLabels=='ev', 'yes'], 0.75)
extremeIcsThreshold = quantile(finalResults[theLabels=='ic', 'no'], 0.75)

evExtremeSamplesId = finalResults[theLabels=='ev' & 
                              finalResults$yes >= extremeEvsThreshold, 'id']
icExtremeSamplesId = finalResults[theLabels=='ic' & 
                              finalResults$no >= extremeIcsThreshold, 'id']

designMatFinal = read.csv('final.design.mat.csv',
                          colClasses = rep('character', 5))

designMatFinal$seq = as.character(DNAStringSet(RNAStringSet(designMatFinal$seq)))

evExtremeSamples = designMatFinal[designMatFinal$id %in% evExtremeSamplesId, ]
icExtremeSamples = designMatFinal[designMatFinal$id %in% icExtremeSamplesId, ]

Preprocess.Write2Fasta(evExtremeSamples, 'ev.extreme.dna.fasta')
Preprocess.Write2Fasta(icExtremeSamples, 'ic.extreme.dna.fasta')

system("../../apps/homer2/cpp/homer2 denovo \
       -i ev.extreme.dna.fasta \
       -b ic.extreme.dna.fasta \
       -len 10 \
       -p 32 \
       -cache 16000 \
       > homerData/homer2.denovo.10.results")
