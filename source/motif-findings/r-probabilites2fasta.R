setwd('/media/pgdrive/sharif/exosomians/')

options(stringsAsFactors = F)
options(scscipen = 999)

library(data.table, quietly = T)
library(Biostrings, quietly = T)

EXOSOMIANS_DIR = '/media/pgdrive/sharif/exosomians'

args = commandArgs(trailingOnly = T)
probsCsv = args[1]
outputName = args[2]
outputType = args[3]

print(outputName)
print(outputType)

# probsCsv = '/media/pgdrive/sharif/exosomians/final.different.purged/final.design.mat.refined.purge80.csv/ExoCNN.ev.extreme.90.refined.probabilities.csv'
# outputName = 'ev.extreme.purged80'
# outputType = 'rna'

probs = fread(probsCsv)

seqs = probs$seq
names(seqs) = probs$id

if(outputType=='rna') {
  seqsFasta = RNAStringSet(seqs)
} else if(outputType == 'dna') {
  seqsFasta = DNAStringSet(RNAStringSet(seqs))
}

outputFasta = sprintf('%s/motif-findings/%s.%s.seqs.fasta', EXOSOMIANS_DIR, outputName, outputType)
writeXStringSet(seqsFasta, outputFasta)
