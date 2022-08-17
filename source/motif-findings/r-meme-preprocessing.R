#### icxAndEvx.purged80 ####

exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'

setwd(exosomiansHomeDir)

library(data.table)
library(Biostrings)

## FASTA files is DNA format
# icAllFastaFile = 'motif-findings/data-motif/ic.all.purged80.dna.seqs.fasta'
icExtremeFastaFile = 'motif-findings/data-motif/ic.extreme.purged80.dna.seqs.fasta'
# evAllFastaFile = 'motif-findings/data-motif/ev.all.purged80.dna.seqs.fasta'
evExtremeFastaFile = 'motif-findings/data-motif/ev.extreme.purged80.dna.seqs.fasta'

# ic = readDNAStringSet(icAllFastaFile)
icx = readDNAStringSet(icExtremeFastaFile)
# ev = readDNAStringSet(evAllFastaFile)
evx = readDNAStringSet(evExtremeFastaFile)

length(icx)
length(evx)

### IC
icxDt = data.table(seq = as.character(icx), id = names(icx))

icxFasta = icxDt$seq
names(icxFasta) = icxDt$id
icxFasta = RNAStringSet(DNAStringSet(icxFasta))

icxFastaFile = 'motif-findings/EXOSOMIANS_DATA/icx.exocnn.purged80.rna.seqs.fasta'
writeXStringSet(icxFasta, icxFastaFile)

### EV
evxDt = data.table(seq = as.character(evx), id = names(evx))
evxFasta = evxDt$seq
names(evxFasta) = evxDt$id
evxFasta = RNAStringSet(DNAStringSet(evxFasta))

evxFastaFile = 'motif-findings/EXOSOMIANS_DATA/evx.exocnn.purged80.rna.seqs.fasta'
writeXStringSet(evxFasta, evxFastaFile)


######################################################## #
#### icAndEv.top.1e4.purged80 ####
### Running on equal number of IC and EV (most promising ones)

exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'

setwd(exosomiansHomeDir)

library(data.table)
library(Biostrings)

allProbsFile = 'final.different.purged/final.design.mat.refined.purge80.csv/ExoCNN.purge80.predictions.csv'

probs = fread(allProbsFile)
icProbs = probs[no>0.5][order(no, decreasing = T)][seq(1e4)]
evProbs = probs[yes>0.5][order(yes, decreasing = T)][seq(1e4)]
head(evProbs)

evProbs[grepl('ev', id), .N]
evProbs[grepl('ev', id), .N, yes>=0.95]
probs[grep('ic', id), .N, no<=0.05]

probs[grepl('ev', id), .N]

### for icAndEv.top.1e4.ex.gt.99.purged80

icProbs = icProbs[, .(id, seq)]
evProbs = evProbs[, .(id, seq)]


### IC
icTop1e4Fasta = icProbs$seq
names(icTop1e4Fasta) = icProbs$id
icTop1e4Fasta = RNAStringSet(icTop1e4Fasta)
icTop1e4FastaFile = 'motif-findings/EXOSOMIANS_DATA/ic.top.1e4.exocnn.purged80.rna.seqs.fasta'
writeXStringSet(icTop1e4Fasta, icTop1e4FastaFile)

### EV
evTop1e4Fasta = evProbs$seq
names(evTop1e4Fasta) = evProbs$id
evTop1e4Fasta = RNAStringSet(evTop1e4Fasta)
evTop1e4FastaFile = 'motif-findings/EXOSOMIANS_DATA/ev.top.1e4.exocnn.purged80.rna.seqs.fasta'
writeXStringSet(evTop1e4Fasta, evTop1e4FastaFile)
