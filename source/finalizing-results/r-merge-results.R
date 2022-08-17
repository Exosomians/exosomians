setwd('/media/pgdrive/sharif/exosomians/final-results')
library(data.table)
library(stringr)
library(Biostrings)

# seqs_mutated_file = 'MutationMaps/top.EVX.greedy.mutation.map.tsv'
# mode = 'top'

.CtrlSeqs2Fasta = function(seqs_mutated_file, mode)
{
  seqs_mutated = fread(seqs_mutated_file)
  seqs_mutated[, id:=str_c('CTRL_', id)]
  
  seqs_mutated_fasta = seqs_mutated$mutated_seq
  names(seqs_mutated_fasta) = seqs_mutated$id
  seqs_mutated_fasta = RNAStringSet(seqs_mutated_fasta)
  seqs_mutated_fasta_file = sprintf('%s/CTRL.fasta', mode)
  
  
  writeXStringSet(seqs_mutated_fasta, seqs_mutated_fasta_file)
}

.CtrlSeqs2Fasta('MutationMaps/top.EVX.greedy.mutation.map.tsv', 'top')
.CtrlSeqs2Fasta('MutationMaps/top.purged.EVX.greedy.mutation.map.tsv', 'top-purged')
.CtrlSeqs2Fasta('MutationMaps/top.stratified.EVX.greedy.mutation.map.tsv', 'top-stratified')
.CtrlSeqs2Fasta('MutationMaps/top.stratified.purged.EVX.greedy.mutation.map.tsv', 'top-stratified-purged')

.MergeResults = function(mode) {
  merge_results_cmd = 'cd %s && cat EVX.fasta CTRL.fasta ICX.fasta REX.fasta RIX.fasta > %s.fasta && cp %s.fasta ..'
  
  system(sprintf(merge_results_cmd, mode, mode, mode))
}

.MergeResults('top')
.MergeResults('top-purged')
.MergeResults('top-stratified')
.MergeResults('top-stratified-purged')

system('zip prioritized-seqs-updated.zip top*')
