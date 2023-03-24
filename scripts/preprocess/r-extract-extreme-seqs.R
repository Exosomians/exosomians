options(stringsAsFactors = F)
options(scipen = 999)

library(Biostrings)

setwd('/media/pgdrive/sharif/exosomians/data/bams/')

designMatShortened = read.csv('final.design.mat.shortened.csv')
evExtremesId = read.csv('deepbind.ev.extreme.predictions.csv')
icExtremesId = read.csv('deepbind.ic.extreme.predictions.csv')

evExtremesDesignMat = designMatShortened[evExtremesId$X+1, ]
icExtremesDesignMat = designMatShortened[icExtremesId$X+1, ]

evExtremesDesignMat$seq = as.character(DNAStringSet(RNAStringSet(evExtremesDesignMat$seq)))
icExtremesDesignMat$seq = as.character(DNAStringSet(RNAStringSet(icExtremesDesignMat$seq)))

set.seed(1)
evExtremesDesignMatSubsampled = evExtremesDesignMat[sample(seq(nrow(evExtremesDesignMat)),
                                                           floor(0.7 * nrow(evExtremesDesignMat))), ]

icExtremesDesignMatSubsampled = icExtremesDesignMat[sample(seq(nrow(icExtremesDesignMat)),
                                                           floor(0.65 * nrow(icExtremesDesignMat))), ]


sum(nchar(icExtremesDesignMatSubsampled$seq))

seqinr::write.fasta(sequences = strsplit(evExtremesDesignMatSubsampled$seq, ''),
                    names = evExtremesDesignMatSubsampled$id,
                    file.out = 'deepbind.ev.extreme.seq.dna.subsampled.fasta')

seqinr::write.fasta(sequences = strsplit(as.character(RNAStringSet(DNAStringSet(evExtremesDesignMatSubsampled$seq))), ''),
                    names = evExtremesDesignMatSubsampled$id,
                    file.out = 'deepbind.ev.extreme.seq.rna.subsampled.fasta')


seqinr::write.fasta(sequences = strsplit(icExtremesDesignMatSubsampled$seq, ''),
                    names = icExtremesDesignMatSubsampled$id,
                    file.out = 'deepbind.ic.extreme.seq.dna.subsampled.fasta')

seqinr::write.fasta(sequences = strsplit(as.character(RNAStringSet(DNAStringSet(icExtremesDesignMatSubsampled$seq))), ''),
                    names = icExtremesDesignMatSubsampled$id,
                    file.out = 'deepbind.ic.extreme.seq.rna.subsampled.fasta')


# evsub = evExtremesDesignMatSubsampled
evsub = evExtremesDesignMat[sample(seq(nrow(evExtremesDesignMat)), 500), ]

icsub = icExtremesDesignMat[sample(seq(nrow(icExtremesDesignMat)), 500), ]

seqinr::write.fasta(sequences = strsplit(evsub$seq, ''),
                    names = evsub$id,
                    file.out = 'deepbind.ev.extreme.seq.dna.subsampled.500.fasta')

seqinr::write.fasta(sequences = strsplit(icsub$seq, ''),
                    names = icsub$id,
                    file.out = 'deepbind.ic.extreme.seq.dna.subsampled.500.fasta')


evxCluste1Ids = read.table('deepbind.ev.extreme.cluster1.ids', header = F)
evxCluste1Ids = evxCluste1Ids$V1

evxCluste1DesignMat = evExtremesDesignMat[evExtremesDesignMat$id%in%evxCluste1Ids, ]
nrow(evxCluste1DesignMat)

seqinr::write.fasta(sequences = strsplit(evxCluste1DesignMat$seq, ''),
                    names = evxCluste1DesignMat$id,
                    file.out = 'deepbind.ev.extreme.cluster1.seq.fasta')


evxCluster2Ids = read.table('deepbind.ev.extreme.cluster2.ids', header = F)
evxCluster2Ids = evxCluster2Ids$V1

evxCluster2DesignMat = evExtremesDesignMat[evExtremesDesignMat$id%in%evxCluster2Ids, ]
nrow(evxCluster2DesignMat)

seqinr::write.fasta(sequences = strsplit(evxCluster2DesignMat$seq, ''),
                    names = evxCluster2DesignMat$id,
                    file.out = 'deepbind.ev.extreme.cluster2.seq.fasta')
