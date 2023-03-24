setwd('/media/pgdrive/sharif/exosomians/')

library(data.table)
library(Biostrings)
library(stringr)

secretionProbs.csvFile = 'final.different.purged/final.design.mat.refined.purgeInf.csv/ExoCNN.purgeInf.seqs.tsv'

secretionProbs = fread(secretionProbs.csvFile)
evProbs = secretionProbs[secretion_prob>0.5][order(secretion_prob, decreasing = T)]

allEvs.fastaFile = 'data/fastas/ev.all.gt0.5.exocnn.rna.fasta'
extremeEvs.fastaFile = 'data/fastas/ev.extreme.gt0.95.exocnn.rna.fasta'

allEvs.fasta = evProbs[, seq]
names(allEvs.fasta) = evProbs[, id]
allEvs.fasta = RNAStringSet(allEvs.fasta)
writeXStringSet(allEvs.fasta, allEvs.fastaFile)

extremeEvs.fasta = evProbs[secretion_prob>0.95 & str_detect(id, 'ev'), seq]
names(extremeEvs.fasta) = evProbs[secretion_prob>0.95 & str_detect(id, 'ev'), id]
extremeEvs.fasta = RNAStringSet(extremeEvs.fasta)
writeXStringSet(extremeEvs.fasta, extremeEvs.fastaFile)
