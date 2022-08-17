EXOSOMIANS_DIR = '/media/pgdrive/sharif/exosomians'
MEME_DIR = '/home/pg/meme/bin' # Worker1 Local

setwd(EXOSOMIANS_DIR)

options(stringsAsFactors = F)
options(scipen = 999)

library(data.table)
library(parallel)
library(Biostrings)
library(stringr)
library(ggplot2)

## Fastas after running r-preprocess on raw bed data
##  BED -> SequenceExtractionAndRefinements -> Write2Fasta
##  Removed seqs containing N
##  Removed long (>50nt) reads)

EV_FASTA = 'data/fastas/ev.final.all.shortened.fasta'
# IC_FASTA = 'data/fastas/ic.final.all.shortened.fasta' ## 3M
IC_FASTA = 'data/fastas/ic.final.all.shortened.downsampled.fasta' ## 3M->160K

## Retrieve previous IC ids
# icProbabilities = fread('predictions/ExoCNN/ExoCNN.final.probabilities.just-ic.csv')
# prevIcIds = icProbabilities$V2
# 
# seqs = readRNAStringSet(IC_FASTA) # read"RNA"StringSet
# seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
# seqs = seqs[id %in% prevIcIds]
# seqsFasta = seqs$seq
# names(seqsFasta) = seqs$id
# seqsFasta = RNAStringSet(seqsFasta)
# writeXStringSet(seqsFasta, 'data/fastas/ic.final.all.shortened.downsampled.fasta')

# evSeqs = Preprocess.RefineSeqs(EV_FASTA, 'ev')
# icSeqs = Preprocess.RefineSeqs(IC_FASTA, 'ic')
# print('Done!')


DUST_CUTOFF = 5
PURGE_SCORE = 'Inf'

GenerateDesignMatrix = function(DUST_CUTOFF, PURGE_SCORE)
{
  LABEL = 'ev'
  FASTA = ifelse(PURGE_SCORE=='Inf',
                 sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.fasta', LABEL, DUST_CUTOFF),
                 sprintf('data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta', LABEL, DUST_CUTOFF, PURGE_SCORE))
  
  if(PURGE_SCORE=='Inf')
  {
    seqs = readRNAStringSet(FASTA)
  } else {
    seqs = readDNAStringSet(FASTA)
    seqs = RNAStringSet(seqs)    
  }
  
  seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
  evSeqs = seqs
  
  LABEL = 'ic'
  FASTA = ifelse(PURGE_SCORE=='Inf',
                 sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.fasta', LABEL, DUST_CUTOFF),
                 sprintf('data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta', LABEL, DUST_CUTOFF, PURGE_SCORE))
  
  if(PURGE_SCORE=='Inf')
  {
    seqs = readRNAStringSet(FASTA)
  } else {
    seqs = readDNAStringSet(FASTA)
    seqs = RNAStringSet(seqs)    
  }

  seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
  icSeqs = seqs
  
  evSeqs$label = 'ev'
  icSeqs$label = 'ic'
  
  allSeqs = rbind(evSeqs, icSeqs)
  
  designMatOutputFile = sprintf('data/csvs/final.design.mat.refined.purge%s.csv', PURGE_SCORE)
  labelMatOutputFile = sprintf('data/csvs/final.label.mat.refined.purge%s.csv', PURGE_SCORE)
  fwrite(allSeqs[, .(id, seq)], file = designMatOutputFile, quote = F)
  fwrite(allSeqs[, .(label)], file = labelMatOutputFile, quote = F)
}


GenerateDesignMatrix(5, 50)
GenerateDesignMatrix(5, 60)
GenerateDesignMatrix(5, 70)
GenerateDesignMatrix(5, 80)
GenerateDesignMatrix(5, 'Inf')

# evSeqs$label = 'ev'
# icSeqs$label = 'ic'
# 
# allSeqs = rbind(evSeqs, icSeqs)
# table(allSeqs$label)
# 
# fwrite(allSeqs[, .(id, seq)], file = 'data/csvs/final.design.mat.refined.csv', quote = F)
# fwrite(allSeqs[, .(label)], file = 'data/csvs/final.label.mat.refined.csv', quote = F)

system('tail data/csvs/final.design.mat.refined.csv')
system('tail data/csvs/final.label.mat.refined.csv')

Preprocess.RefineSeqs = function(fastaFile, label)
{
  
  ## STAGE 1: Remove duplicated seqs
  ## STAGE 2: Remove substrings
  ## STAGE 3: Remove/Mask low complex seqs
  ## STAGE 4: Remove highly similar seqs
  
  FASTA = fastaFile
  # FASTA = IC_FASTA
  LABEL = label
  # LABEL = 'ic'
  
  seqs = readRNAStringSet(FASTA) # read"RNA"StringSet
  seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
  
  # nrow(seqs)
  ## 40221
  
  setindex(seqs, seq)
  #### STAGE 1: Remove duplicated seqs ####
  seqs = unique(seqs, by = 'seq')
  
  # nrow(seqs)
  ## 39080
  
  #### STAGE 2: Remove substrings ####
  seqsArray = as.character(seqs$seq)
  names(seqsArray) = seqs$id
  
  seqsList = as.list(seqs$seq)
  
  ## Detect substrings
  isSubstring = mclapply(seqsList, function(aSeq)
    ifelse(sum(str_detect(seqsArray, aSeq), na.rm = T) > 1, TRUE, FALSE),
    mc.cores = detectCores() - 2)
  
  isSubstring = unlist(isSubstring)
  # table(isSubstring)
  
  ## Remove substrings
  seqs = seqs[!isSubstring]
  
  # nrow(seqs)
  ## 38144
  
  #### STAGE 3: Remove/Mask low complex seqs ####
  
  ## We should do it using MEME dust
  
  seqsFasta = seqs$seq
  names(seqsFasta) = seqs$id
  seqsFasta = RNAStringSet(seqsFasta)
  
  writeXStringSet(seqsFasta, filepath = sprintf('data/fastas/%s.preprocessed.stage2.fasta', LABEL))
  
  ## BASH
  # dustCutoffs = c(seq(5), 10, 15, 20)

  
  ##### TO BE uncommented!
  dustCutoffs = 5
  seqsStage3 = lapply(dustCutoffs, function(aCutoff)
  {
    DUST_CUTOFF = aCutoff
    system(sprintf('%s/dust data/fastas/%s.preprocessed.stage2.fasta %s > data/fastas/%s.preprocessed.stage3.dust%s.fasta',
                   MEME_DIR,
                   LABEL,
                   DUST_CUTOFF,
                   LABEL,
                   DUST_CUTOFF))

    seqs = readRNAStringSet(sprintf('data/fastas/%s.preprocessed.stage3.dust%s.fasta', LABEL, DUST_CUTOFF)) # read"RNA"StringSet
    seqs = data.table(id = names(seqs), seq = as.character(seqs))
    setnames(seqs, 'seq', paste0('seq_dust', DUST_CUTOFF))
    seqs
  })

  names(seqsStage3) = paste0('dust', dustCutoffs)
  
  # nrow(seqs)
  ## 38144
  ## No seqs removed, just masked
  
  
  ### Let's compare and test different dust cut offs
  
  # seqsStage2 = readRNAStringSet(sprintf('data/fastas/%s.preprocessed.stage2.fasta', LABEL)) 
  # seqsStage2 = data.table(id = names(seqsStage2), seq_normal = as.character(seqsStage2))
  # 
  # ## seqsStage2, a data.table
  # ## seqsStage3, a list of data.tables
  # 
  # seqsStage2And3 = c(list(seqsStage2), seqsStage3)
  # seqsComparison = Reduce(function(x, y) merge(x, y, by = 'id', all.x = T), seqsStage2And3)
  # 
  # seqsStage3Cols = paste0('seq_', names(seqsStage3))
  # seqsComparison[, (seqsStage3Cols):= lapply(.SD, function(seq_dustX)
  #   {
  #   seq_dustX[seq_dustX == seqsComparison[, seq_normal]] = NA
  #   seq_dustX
  #   }), .SDcols = seqsStage3Cols]
  # 
  # 
  # View(seqsComparison)
  # seqsCounts = apply(seqsComparison, 2, function(aCol) length(aCol[is.na(aCol)]))
  # seqsCounts = data.table(dust = names(seqsCounts[-seq(2)]), count = seqsCounts[-seq(2)])
  # 
  # 
  # ggplot(seqsCounts) +
  #   aes(x = reorder(dust, count), y = count, group = 1) +
  #   geom_line(linetype = 'dashed') +
  #   geom_point()
  
  
  ### Selected DUST_CUTOFFs -> 5, 10
  DUST_CUTOFF = 5
  
  seqs = readRNAStringSet(sprintf('data/fastas/%s.preprocessed.stage3.dust%s.fasta', LABEL, DUST_CUTOFF)) # read"RNA"StringSet
  seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
  
  ### Remove seqs that have been masked by dust (have N)
  seqs = seqs[!str_detect(seqs$seq, 'N')]

  # nrow(seqs)
  ## 33082
  
  #### STAGE 4: Remove highly similar seqs ####
  seqsFasta = seqs$seq
  names(seqsFasta) = seqs$id
  seqsFasta = RNAStringSet(seqsFasta)
  
  writeXStringSet(seqsFasta, filepath = sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.fasta', LABEL, DUST_CUTOFF))
  
  
  seqsFasta = DNAStringSet(seqsFasta)
  writeXStringSet(seqsFasta, filepath = sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.dna.fasta', LABEL, DUST_CUTOFF))
  
  # purgeScores = c(1, 5, seq(10, 100, 10))
  # purgeScores = c(50, 60, 70, 80) ## 50 is already done!
  print('Purge started!')
  purgeScores = c(60, 70, 80) 
  seqsStage4 = mclapply(purgeScores, function(aScore)
  {
    PURGE_SCORE = aScore
    system(sprintf('%s/purge data/fastas/%s.preprocessed.stage3.dust%s.final.dna.fasta %s -n -o > data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta',
                   MEME_DIR,
                   LABEL,
                   DUST_CUTOFF,
                   PURGE_SCORE,
                   LABEL,
                   DUST_CUTOFF,
                   PURGE_SCORE))
    
    seqs = readDNAStringSet(sprintf('data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta', LABEL, DUST_CUTOFF, PURGE_SCORE))
    seqs = RNAStringSet(seqs)
    seqs = data.table(id = names(seqs), seq = as.character(seqs))
    setnames(seqs, 'seq', paste0('seq_dust', DUST_CUTOFF, '_purge', PURGE_SCORE))
    seqs
  }, mc.cores = detectCores() - 2)
  
  names(seqsStage4) = paste0('dust', DUST_CUTOFF, '_purge', purgeScores)
  
  ### Let's compare and test different purge scores
  
  # seqsStage3 = readRNAStringSet(sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.fasta', LABEL, DUST_CUTOFF)) 
  # seqsStage3 = data.table(id = names(seqsStage3), seq_normal = as.character(seqsStage3))
  # 
  # ## seqsStage3, a data.table
  # ## seqsStage4, a list of data.tables
  # 
  # seqsStage3And4 = c(list(seqsStage3), seqsStage4)
  # seqsComparison = Reduce(function(x, y) merge(x, y, by = 'id', all.x = T), seqsStage3And4)
  # 
  # View(seqsComparison)
  # seqsCounts = apply(seqsComparison, 2, function(aCol) length(aCol[!is.na(aCol)]))
  # seqsCounts = data.table(purge = names(seqsCounts[-seq(2)]), count = seqsCounts[-seq(2)])
  # 
  # ggplot(seqsCounts) +
  #   aes(x = reorder(purge, count), y = count, group = 1) +
  #   geom_line(linetype = 'dashed') +
  #   geom_point()
  
  
  ### Purge Tuning
  ## Using MEME motif finder
  
  # seqsNoPurge = fread('data/fastas/purge-tuning/meme-no-purge-no-header.ids', header = F)
  # seqsNoPurge = seqsNoPurge[, c(2, 3)]
  # colnames(seqsNoPurge) = c('id', 'seq')
  # seqsNoPurgeFasta = seqsNoPurge$seq
  # names(seqsNoPurgeFasta) = seqsNoPurge$id
  # seqsNoPurgeFasta = RNAStringSet(seqsNoPurgeFasta)
  # 
  # writeXStringSet(seqsNoPurgeFasta, 'data/fastas/purge-tuning/meme-no-purge-no-header.fasta')
  
  ### Before, DUST_CUTOFF -> 5
  ### Selected PURGE_SCORE -> 50
  PURGE_SCORE = 50
  
  seqs = readDNAStringSet(sprintf('data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta', LABEL, DUST_CUTOFF, PURGE_SCORE))
  seqs = RNAStringSet(seqs)
  seqs = data.table(id = names(seqs), seq = as.character(seqs), width = width(seqs))
  
  # nrow(seqs)
  ## 5757
  
  seqs
}



### Let's compare and test different purge scores

LABEL = 'ev'
DUST_CUTOFF = 5

seqsStage3 = readRNAStringSet(sprintf('data/fastas/%s.preprocessed.stage3.dust%s.final.fasta', LABEL, DUST_CUTOFF))
seqsStage3 = data.table(id = names(seqsStage3), seq_normal = as.character(seqsStage3))

## seqsStage3, a data.table
## seqsStage4, a list of data.tables

purgeScores = rev(c(seq(5), seq(10, 100, 10)))
seqsStage4 = mclapply(purgeScores, function(aScore)
{
  PURGE_SCORE = aScore
  seqs = readDNAStringSet(sprintf('data/fastas/%s.preprocessed.stage4.dust%s.purge%s.fasta', LABEL, DUST_CUTOFF, PURGE_SCORE))
  seqs = RNAStringSet(seqs)
  seqs = data.table(id = names(seqs), seq = as.character(seqs))
  setnames(seqs, 'seq', paste0('seq_dust', DUST_CUTOFF, '_purge', PURGE_SCORE))
  seqs
}, mc.cores = detectCores() - 2)

seqsStage3And4 = c(list(seqsStage3), seqsStage4)
seqsComparison = Reduce(function(x, y) merge(x, y, by = 'id', all.x = T), seqsStage3And4)

View(seqsComparison)
seqsCounts = apply(seqsComparison, 2, function(aCol) length(aCol[!is.na(aCol)]))
seqsCounts = data.table(purge = names(seqsCounts[-seq(1)]), count = seqsCounts[-seq(1)])

fwrite(seqsComparison, 'data/csvs/ev.extreme.refined.purged.seqs-comparison.csv')
fwrite(seqsCounts, 'data/csvs/ev.extreme.refined.purged.seqs-counts.csv')

ggplot(seqsCounts) +
  aes(x = reorder(purge, count), y = count, group = 1) +
  geom_line(linetype = 'dashed') +
  geom_point()
