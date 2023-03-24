exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'
fireHomeDir = 'apps/FIRE-1.1a'

setwd(exosomiansHomeDir)

library(data.table)
library(Biostrings)


## FASTA files is DNA format
icAllFastaFile = 'motif-findings/data-motif/ic.all.purged80.dna.seqs.fasta'
icExtremeFastaFile = 'motif-findings/data-motif/ic.extreme.purged80.dna.seqs.fasta'
evAllFastaFile = 'motif-findings/data-motif/ev.all.purged80.dna.seqs.fasta'
evExtremeFastaFile = 'motif-findings/data-motif/ev.extreme.purged80.dna.seqs.fasta'

ic = readDNAStringSet(icAllFastaFile)
icx = readDNAStringSet(icExtremeFastaFile)
ev = readDNAStringSet(evAllFastaFile)
evx = readDNAStringSet(evExtremeFastaFile)

length(ic)
length(ev)

icDt = data.table(seq = as.character(ic), id = names(ic))
evDt = data.table(seq = as.character(ev), id = names(ev))

icDt[, loc:=1]
evDt[, loc:=2]

icDt[id %in% names(icx), loc:=0]
evDt[id %in% names(evx), loc:=3]

icAndEvDt = rbind(icDt, evDt)

icAndEvFasta = icAndEvDt$seq
names(icAndEvFasta) = icAndEvDt$id

icAndEvFasta = DNAStringSet(icAndEvFasta)
icAndEvFastaFile = 'motif-findings/data-motif/icAndEv.all.purged80.dna.seqs.fasta'
writeXStringSet(icAndEvFasta, icAndEvFastaFile)

system(sprintf('head %s', icAndEvFastaFile))

icAndEvExp = icAndEvDt[, .(id, loc)][order(loc)]
setnames(icAndEvExp, 'id', 'smRNA_id')

icAndEvExpFile = 'motif-findings/data-motif/icAndEv.all.purged80.expfile.txt'
fwrite(icAndEvExp, icAndEvExpFile, sep = '\t')

system(sprintf('head %s', icAndEvExpFile))

args = list()
args$exptype = 'discrete'

fireCommand = sprintf('export FIREDIR=%s && perl $FIREDIR/fire.pl --expfile=%s  --species=human --exptype=%s --fastafile_rna=%s --oribiasonly=0 --shuffle=1000 --jn_t=0 --dodna=0 --dodnarna=0 --nodups=1 --shuffle_mifind=1000 --shuffle_midist=1000',
                      paste(exosomiansHomeDir, fireHomeDir, sep = '/'),
                      icAndEvExpFile,
                      args$exptype,
                      icAndEvFastaFile)

system(fireCommand) # Running out of memory with these numbers of sequences
# 127120 ic + 28213 ev = 155333 seqs


## ERROR! No data file for human.

# icAndEvFastaFile = paste(fireHomeDir, 'icAndEv.all.purged80.dna.seqs.fasta', sep = '/')
# writeXStringSet(icAndEvFasta, icAndEvFastaFile)
# 
# icAndEvExpFile = paste(fireHomeDir, 'icAndEv.all.purged80.expfile.txt', sep = '/')
# fwrite(icAndEvExp, icAndEvExpFile, sep = '\t')
# 
# fireCommand = sprintf('export FIREDIR=%s && perl $FIREDIR/fire.pl --expfile=%s  --species=human --exptype=%s --fastafile_rna=%s --oribiasonly=0 --shuffle=1000 --jn_t=0 --dodna=0 --dodnarna=0 --nodups=1 --shuffle_mifind=1000 --shuffle_midist=1000',
#                       paste(exosomiansHomeDir, fireHomeDir, sep = '/'),
#                       gsub('.*/', '', icAndEvExpFile),
#                       args$exptype,
#                       gsub('.*/', '', icAndEvFastaFile))
# 
# system(fireCommand)

######################################################## #
#### icAndEv.top.1e4.purged80 and icAndEv.top.1e4.ex.gt.99.purged80 ####
### Running on equal number of IC and EV (most promising ones)


exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'
fireHomeDir = 'apps/FIRE-1.1a'

setwd(exosomiansHomeDir)

library(data.table)
library(Biostrings)

allProbsFile = 'final.different.purged/final.design.mat.refined.purge80.csv/ExoCNN.purge80.predictions.csv'

probs = fread(allProbsFile)
icProbs = probs[no>0.5][order(no, decreasing = T)][seq(1e4)]
evProbs = probs[yes>0.5][order(yes, decreasing = T)][seq(1e4)]

### for icAndEv.top.1e4.ex.gt.99.purged80
evProbs[yes>=quantile(yes, 0.9), summary(yes)]
icProbs[no>=0.9999, .N]

icProbs[, loc:=1]
evProbs[, loc:=2]
icProbs[no>=0.9999, loc:=0]
evProbs[yes>=0.99, loc:=3]

icProbs = icProbs[, .(id, seq, loc)]
evProbs = evProbs[, .(id, seq, loc)]

icAndEvProbs = rbind(icProbs, evProbs)
icAndEvProbs[, .N, loc][order(loc)]

icAndEvFasta = icAndEvProbs$seq
names(icAndEvFasta) = icAndEvProbs$id
icAndEvFasta = DNAStringSet(RNAStringSet(icAndEvFasta))
icAndEvFastaFile = 'motif-findings/data-motif/icAndEv.top.1e4.ex.gt.99.purged80.dna.seqs.fasta'
writeXStringSet(icAndEvFasta, icAndEvFastaFile)

icAndEvExp = icAndEvProbs[, .(id, loc)][order(loc)]
setnames(icAndEvExp, 'id', 'smRNA_id')

icAndEvExpFile = 'motif-findings/data-motif/icAndEv.top.1e4.ex.gt.99.purged80.expfile.txt'
fwrite(icAndEvExp, icAndEvExpFile, sep = '\t')

args = list()
args$exptype = 'discrete'

fireCommand = sprintf('export FIREDIR=%s && perl $FIREDIR/fire.pl --expfile=%s  --species=human --exptype=%s --fastafile_rna=%s --oribiasonly=0 --shuffle=1000 --jn_t=0 --dodna=0 --dodnarna=0 --nodups=1 --shuffle_mifind=1000 --shuffle_midist=1000',
                      paste(exosomiansHomeDir, fireHomeDir, sep = '/'),
                      icAndEvExpFile,
                      args$exptype,
                      icAndEvFastaFile)

# system(fireCommand)


################################################# #
#### icxAndEvx.all.purged80 ####
### Running just on extreme sequences

exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'
fireHomeDir = 'apps/FIRE-1.1a'

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

icxDt = data.table(seq = as.character(icx), id = names(icx))
evxDt = data.table(seq = as.character(evx), id = names(evx))

icxDt[, loc:=0]
evxDt[, loc:=3]

icxAndEvxDt = rbind(icxDt, evxDt)

icxAndEvxFasta = icxAndEvxDt$seq
names(icxAndEvxFasta) = icxAndEvxDt$id

icxAndEvxFasta = DNAStringSet(icxAndEvxFasta)
icxAndEvxFastaFile = 'motif-findings/data-motif/icxAndEvx.all.purged80.dna.seqs.fasta'
writeXStringSet(icxAndEvxFasta, icxAndEvxFastaFile)

system(sprintf('head %s', icxAndEvxFastaFile))

icxAndEvxExp = icxAndEvxDt[, .(id, loc)][order(loc)]
setnames(icxAndEvxExp, 'id', 'smRNA_id')

icxAndEvxExpFile = 'motif-findings/data-motif/icxAndEvx.all.purged80.expfile.txt'
fwrite(icxAndEvxExp, icxAndEvxExpFile, sep = '\t')

system(sprintf('head %s', icxAndEvxExpFile))

args = list()
args$exptype = 'discrete'

fireCommand = sprintf('export FIREDIR=%s && perl $FIREDIR/fire.pl --expfile=%s  --species=human --exptype=%s --fastafile_rna=%s --oribiasonly=0 --shuffle=1000 --jn_t=0 --dodna=0 --dodnarna=0 --nodups=1 --shuffle_mifind=1000 --shuffle_midist=1000',
                      paste(exosomiansHomeDir, fireHomeDir, sep = '/'),
                      icxAndEvxExpFile,
                      args$exptype,
                      icxAndEvxFastaFile)

system(fireCommand) # Running out of memory with these numbers of sequences

#################################################### #
### Running all all sequences after purge (remove highly similar seqs)

exosomiansHomeDir = '/media/pgdrive/sharif/exosomians'
fireHomeDir = 'apps/FIRE-1.1a'

setwd(exosomiansHomeDir)

library(data.table)
library(Biostrings)

PURGE_SCORE = 70
outputPrefix = 'icAndEv.all.purged80'

## FASTA files is DNA format and is after being PURGED!
icAllFastaFile = paste0('motif-findings/data-motif/ic.all.purged80.dna.seqs.fasta.b', PURGE_SCORE)
evAllFastaFile = paste0('motif-findings/data-motif/ev.all.purged80.dna.seqs.fasta.b', PURGE_SCORE)

allProbsFile = 'final.different.purged/final.design.mat.refined.purge80.csv/ExoCNN.purge80.predictions.csv'

ic = readDNAStringSet(icAllFastaFile)
ev = readDNAStringSet(evAllFastaFile)
probs = fread(allProbsFile)

icDt = data.table(seq = as.character(ic), id = names(ic))
evDt = data.table(seq = as.character(ev), id = names(ev))

icProbs = probs[icDt, on = c('id')][no>=0.9]
evProbs = probs[evDt, on = c('id')][yes>=0.9]


evProbs[, summary(yes)]
evProbs[yes>=0.9, .N]
evProbs[yes>=0.99, .N]
icProbs[no>=0.9, .N]
icProbs[no>=0.99, .N]


### for numOfClasses==4
# icProbs[, loc:=1]
# evProbs[, loc:=2]
# 
# icProbs[no>=0.99, loc:=0]
# evProbs[yes>=0.99, loc:=3]

### for numOfClasess==2
icProbs[, loc:=0]
evProbs[, loc:=3]


icProbs = icProbs[, .(id, seq, loc)]
evProbs = evProbs[, .(id, seq, loc)]

icAndEvProbs = rbind(icProbs, evProbs)
icAndEvProbs[, .N, loc][order(loc)]
numOfClasses = icAndEvProbs[, length(unique(loc))]

### FASTA
icAndEvFasta = icAndEvProbs$seq
names(icAndEvFasta) = icAndEvProbs$id

icAndEvFasta = DNAStringSet(RNAStringSet(icAndEvFasta))

icAndEvFastaFile = sprintf('motif-findings/data-motif/%s.dna.seqs.b%s.c%s.fasta',
                           outputPrefix,
                           PURGE_SCORE,
                           numOfClasses)

writeXStringSet(icAndEvFasta, icAndEvFastaFile)

### EXPFILE
icAndEvExp = icAndEvProbs[, .(id, loc)][order(loc)]
setnames(icAndEvExp, 'id', 'smRNA_id')

icAndEvExpFile = sprintf('motif-findings/data-motif/%s.expfile.b%s.c%s.txt',
                         outputPrefix,
                         PURGE_SCORE,
                         numOfClasses)

fwrite(icAndEvExp, icAndEvExpFile, sep = '\t')

### FIRE
targetDir = sprintf('%s/EXOSOMIANS_DATA/%s.b%s.c%s',
                    fireHomeDir,
                    outputPrefix,
                    PURGE_SCORE,
                    numOfClasses)

dir.create(targetDir)

file.copy(from = icAndEvFastaFile,
          to = sprintf('%s/%s',
                       targetDir,
                       gsub('.*/', '', icAndEvFastaFile)))

file.copy(from = icAndEvExpFile,
          to = sprintf('%s/%s',
                       targetDir,
                       gsub('.*/', '', icAndEvExpFile)))


fireCommand = sprintf('export FIREDIR=%s && cd %s && perl $FIREDIR/fire.pl --expfile=%s  --species=human 1--exptype=discrete --fastafile_rna=%s --oribiasonly=0 --shuffle=1000 --jn_t=0 --dodna=0 --dodnarna=0 --nodups=1 --shuffle_mifind=1000 --shuffle_midist=1000',
                      paste(exosomiansHomeDir, fireHomeDir, sep = '/'),
                      targetDir,
                      gsub('.*/', '', icAndEvExpFile),
                      gsub('.*/', '', icAndEvFastaFile))

fireCommand
# system(fireCommand)
