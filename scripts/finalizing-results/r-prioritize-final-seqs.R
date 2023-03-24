EXOSOMIANS_DIR = '/media/pgdrive/sharif/exosomians'
N = 2000/5 # Number of candidate seqs
PURGE_DIR = '~/meme/bin/purge'
PURGE_SCORE = 50
  
# dir.create(sprintf('%s/final-results/top', EXOSOMIANS_DIR), recursive = T)
# dir.create(sprintf('%s/final-results/top-stratified', EXOSOMIANS_DIR))
# dir.create(sprintf('%s/final-results/top-purged', EXOSOMIANS_DIR))
# dir.create(sprintf('%s/final-results/top-stratified-purged', EXOSOMIANS_DIR))

setwd(EXOSOMIANS_DIR)

library(data.table)
library(stringr)
library(Biostrings)

PrioritizeResults = function(input_seqs_file, label, is_random = F)
{
  input_seqs = fread(input_seqs_file)
  
  head(input_seqs)
  summary(nchar(input_seqs$seq))
  summary(input_seqs$secretion_prob)
  
  #### Extract seqs with highest secretion probabilities ####
  
  top_seqs = input_seqs[order(-secretion_prob)][seq(N*10)]
  head(top_seqs)
  
  chr_len = fread('scripts/hg38-chr-lenghts.txt', col.names = c('chr', 'len'))
  chr_len[, frac:=len/sum(len)]
  chr_len[, num:=ceiling((N*10)*frac)]
  
  # sum(chr_len$num)
  
  #### Preparing for the purge tool ####
  ### top_seqs
  
  .Seq2Fasta = function(sequences)
  {
    fasta_seqs = sequences$seq
    names(fasta_seqs) = sequences$id
    fasta_seqs = RNAStringSet(fasta_seqs)
    fasta_seqs = DNAStringSet(fasta_seqs)
    fasta_seqs
  }
  
  top_seqs_fasta = .Seq2Fasta(top_seqs)
  top_seqs_fasta_file = str_replace(input_seqs_file, '.tsv', '.prioritized.dna.fasta')
  writeXStringSet(top_seqs_fasta, top_seqs_fasta_file)
  top_seqs_purged_fasta_file = str_replace(top_seqs_fasta_file, '.fasta', '.purged.fasta')
  
  #### Applying purge to remove highly similar seqs ####
  
  purge_cmd = sprintf('%s %s %s -n -o > %s',
                      PURGE_DIR,
                      top_seqs_fasta_file,
                      PURGE_SCORE,
                      top_seqs_purged_fasta_file)
  
  system(purge_cmd)
  
  
  #### Make the number of seqs equal to N ####
  
  top_seqs_final = top_seqs[seq(N)]
  
  .Fasta2Seq = function(fasta_file)
  {
    seqs_fasta = readDNAStringSet(fasta_file)
    seqs_fasta = RNAStringSet(seqs_fasta)
    seqs = data.table(id = names(seqs_fasta), seq = as.character(seqs_fasta))
    seqs = input_seqs[seqs, on = .(id, seq)]
    seqs
  }
  
  top_seqs_purged = .Fasta2Seq(top_seqs_purged_fasta_file)
  top_seqs_purged_final = top_seqs_purged[order(-secretion_prob)][seq(N)]
  
  #### Final results ####
  ### top_seqs
  ### top_seqs_purged
  
  .FinalizeResults = function(seqs, mode, label) {
    if(!is_random)
    {
      genome_coords = as.data.table(str_split(seqs$id, pattern = '_', simplify = T)[, seq(4)])
      names(genome_coords) = c('chr', 'start', 'end', 'strand')
      seqs = cbind(seqs, genome_coords)
    }
    seqs[, id:=str_c(label, '_', id)]
    
    ### Writing Fasta
    seqs_fasta = seqs$seq
    names(seqs_fasta) = seqs$id
    seqs_fasta = RNAStringSet(seqs_fasta)
    writeXStringSet(seqs_fasta, filepath = sprintf('final-results/%s/%s.fasta', mode, label))
    
    ### Writing CSV
    fwrite(seqs, sprintf('final-results/%s/%s.csv', mode, label), col.names = F)
  }
  
  .FinalizeResults(top_seqs_final, 'top', label)
  .FinalizeResults(top_seqs_purged_final, 'top-purged', label)
  
  
  #### Stratified results (Just for original (non-randomly) generated seqs)
  
  if(!is_random)
  {
    #### Extract seqs with highest prob. stratified by chr (with respect to chr length) ####
    
    input_seqs[, chr:=str_split(id, pattern = '_', n = 2, simplify = T)[, 1]]
    # input_seqs[, .N, chr]
    
    chr_len[, possible_num:=pmin(num, input_seqs[, .N, chr][, N])]
    # sum(chr_len$possible_num) # It ok!
    
    .GetNumOfSeqs = function(chromosome) {
      # browser()
      chromosome = str_remove(chromosome, 'chr')
      chr_len[chr==chromosome, possible_num]
    }
    
    top_seqs_stratified = input_seqs[, .SD[seq(.GetNumOfSeqs(chr))], chr]
    
    #### Preparing for the purge tool ####
    ### top_seqs_stratified
    
    top_seqs_stratified_fasta = .Seq2Fasta(top_seqs_stratified)
    top_seqs_stratified_fasta_file = str_replace(input_seqs_file, '.tsv', '.prioritized.stratified.dna.fasta')
    writeXStringSet(top_seqs_stratified_fasta, top_seqs_stratified_fasta_file)
    top_seqs_stratified_purged_fasta_file = str_replace(top_seqs_stratified_fasta_file, '.fasta', '.purged.fasta')
    
    #### Applying purge to remove highly similar seqs ####
    
    purge_cmd = sprintf('%s %s %s -n -o > %s',
                        PURGE_DIR,
                        top_seqs_stratified_fasta_file,
                        PURGE_SCORE,
                        top_seqs_stratified_purged_fasta_file)
    
    system(purge_cmd)
    
    #### Make the number of seqs equal to N ####
    
    top_seqs_stratified_final = top_seqs_stratified[, .SD[seq(ceiling(.GetNumOfSeqs(chr)/10))], chr][sample(seq(.N), N)]
    top_seqs_stratified_final[, chr:=NULL]
    
    top_seqs_stratified_purged = .Fasta2Seq(top_seqs_stratified_fasta_file)
    top_seqs_stratified_purged_final = top_seqs_stratified_purged[, .SD[seq(ceiling(.GetNumOfSeqs(chr)/10))], chr][sample(seq(.N), N)]
    top_seqs_stratified_purged_final[, chr:=NULL]
    
    #### Final results ####
    ### top_seqs_stratified
    ### top_seqs_stratified_purged
    
    .FinalizeResults(top_seqs_stratified_final, 'top-stratified', label)
    .FinalizeResults(top_seqs_stratified_purged_final, 'top-stratified-purged', label)  
  }
}


## EVX
PrioritizeResults('final.different.purged/final.design.mat.refined.purgeInf.csv/ExoCNN.purgeInf.ev.extreme.95.refined.seqs.tsv', 'EVX')

## ICX
PrioritizeResults('final.different.purged/final.design.mat.refined.purgeInf.csv/ExoCNN.purgeInf.ic.extreme.95.refined.seqs.tsv', 'ICX')

## REX
PrioritizeResults('final.different.purged/final.design.mat.refined.purgeInf.csv/RandomSequences/ExoCNN.ev.extreme.95.random.sequences.3.tsv', 'REX', is_random = T)

## RIX
PrioritizeResults('final.different.purged/final.design.mat.refined.purgeInf.csv/RandomSequences/ExoCNN.ic.extreme.95.random.sequences.3.tsv', 'RIX', is_random = T)
