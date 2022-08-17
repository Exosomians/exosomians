setwd('/media/pgdrive/sharif/exosomians/data')

options(stringsAsFactors = F)
source('scripts/r-functions.R')

library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)

hg38 = BSgenome.Hsapiens.UCSC.hg38

#### Extracting features using annotation files ####
EV_LABEL = 'ev'
EV_BED = 'beds/ev.final.bed'
EV_FASTA = 'fastas/ev.final.fasta'

IC_LABEL = 'ic'
IC_BED = 'beds/ic.final.bed'
IC_FASTA = 'fastas/ic.final.fasta'


#### EV ####
# EV_LABEL = 'ev'
# EV_BED = 'beds/ev.final.in.fimo.merged.0.9.uniq.bed'
# ev.seqs-related-to-encode-srsf1
evDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = EV_BED,
                                                          label = EV_LABEL,
                                                          genome = hg38)
head(evDesignMat)
write.csv(evDesignMat[, c('id', 'seq')], quote = F, row.names = F, file = 'csvs/ev.seqs-related-to-encode-srsf1.csv')

# MAX_SEQ_LEN = 50
# evDesignMatFilteredByLength = evDesignMat[evDesignMat$length<=MAX_SEQ_LEN, ]
# EV_SHORTENED_FASTA = 'fastas/ev.final.all.shortened.fasta'
# Preprocess.Write2Fasta(designMat = evDesignMatFilteredByLength,
#                        outputFile = EV_SHORTENED_FASTA)

EV_FASTA = 'fastas/ev.final.all.fasta'
# EV_FASTA = 'fastas/ev.seqs-related-to-encode-srsf1.fasta'

Preprocess.Write2Fasta(designMat = evDesignMat,
                       outputFile = EV_FASTA)

system("bash -xv main-vienna.sh 'ev.final.fasta'")

evDesignMat = Preprocess.AddDotBracketSeq(evDesignMat,
                                          dotBracketFasta = gsub('.fasta',
                                                                 '.dotbracket.fasta',
                                                                 EV_FASTA))
write.csv(evDesignMat, file = paste0(EV_LABEL, '_DesignMatrix_SeqPlusDotBracket.csv'),
          quote = F, row.names = F)

# evDesignMat = read.csv('ev_DesignMatrix_SeqPlusDotBracket_forgi.csv')
# head(evDesignMat)

#### IC ####

# icDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = IC_BED,
#                                                         label = IC_LABEL,
#                                                         genome = hg38,
#                                                         doDownsample = T,
#                                                         numOfSamples = nrow(evDesignMat),
#                                                         seed = 1)

# IC_LABEL = 'ic'
# IC_BED = 'beds/ic.final.in.fimo.merged.uniq.plus.predictions.bed'
# ic.seqs-related-to-encode-srsf1
icDesignMat = Preprocess.SequenceExtractionAndRefinements(bedFilePath = IC_BED,
                                                          label = IC_LABEL,
                                                          genome = hg38)
quantile(icDesignMat$length, seq(0.95, 1, 0.005))

# MAX_SEQ_LEN = 50
# icDesignMatFilteredByLength = icDesignMat[icDesignMat$length<=MAX_SEQ_LEN, ]
# IC_SHORTENED_FASTA = 'fastas/ic.final.all.shortened.fasta'
# Preprocess.Write2Fasta(designMat = icDesignMatFilteredByLength,
#                        outputFile = IC_SHORTENED_FASTA)

# designMatFilteredByLength = rbind(evDesignMatFilteredByLength[, c('id', 'seq')],
#                                   icDesignMatFilteredByLength[, c('id', 'seq')])
# labelMatFilteredByLength = data.frame(label = c(evDesignMatFilteredByLength[, 'label'],
#                                                icDesignMatFilteredByLength[, 'label']))

write.csv(designMatFilteredByLength, file = 'csvs/final.design.mat.all.shortened.csv',
          quote = F, row.names = F, col.names = T)
write.csv(labelMatFilteredByLength, file = 'csvs/final.label.mat.all.shortened.csv',
          quote = F, row.names = F, col.names = T)


# icDesignMat = Preprocess.Downsample(icDesignMat, evDesignMat, seed = 1)

IC_FASTA = 'fastas/ic.final.all.fasta'
# IC_FASTA = 'fastas/ic.seqs-related-to-encode-srsf1.fasta'

Preprocess.Write2Fasta(designMat = icDesignMat,
                       outputFile = IC_FASTA)

system("bash -xv main-vienna.sh 'ic.final.fasta'")

icDesignMat = Preprocess.AddDotBracketSeq(icDesignMat,
                                          dotBracketFasta = gsub('.fasta',
                                                                 '.dotbracket.fasta',
                                                                 IC_FASTA))
write.csv(icDesignMat, file = paste0(IC_LABEL, '_DesignMatrix_SeqPlusDotBracket.csv'),
          quote = F, row.names = F)

# system("bash -xv main-forgi.sh ic_DesignMatrix_SeqPlusDotBracket.csv")

options(scipen=999)

evDesignMatFinal = read.csv('ev_DesignMatrix_SeqPlusDotBracket_forgi.csv',
                            colClasses = rep('character', 9))
icDesignMatFinal = read.csv('ic_DesignMatrix_SeqPlusDotBracket_forgi.csv',
                            colClasses = rep('character', 9))

evDesignMatFinal$length = as.numeric(evDesignMatFinal$length)
icDesignMatFinal$length = as.numeric(icDesignMatFinal$length)

designMatFinal = rbind(evDesignMatFinal, icDesignMatFinal)

theLabels =  designMatFinal$label
theLabels[theLabels=='ev'] = 'YES'
theLabels[theLabels=='ic'] = 'NO'
theLabels = data.frame(label = theLabels)

designMatFinal = designMatFinal[, c('id', 'seq', 'dotbracket',
                                    'element_string', 'element_string_number')]

write.csv(designMatFinal, file = 'final.design.mat.csv', quote = F, row.names = F)
write.csv(theLabels, file = 'final.label.mat.csv', quote = F, row.names = F)
