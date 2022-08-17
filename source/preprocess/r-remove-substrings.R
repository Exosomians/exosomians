setwd('/media/pgdrive/sharif/exosomians/exosomians')

options(stringsAsFactors = F)
options(scipen = 999)

library(data.table)
library(stringr)
library(parallel)

evPredictionsInputFile = 'predictions/ExoCNN/ExoCNN.ev.extreme.90.probabilities.csv'
icPredictionsInputFile = 'predictions/ExoCNN/ExoCNN.ic.extreme.90.probabilities.csv'
evPredictionsOutputFile = 'predictions/ExoCNN/ExoCNN.ev.extreme.90.unique.probabilities.csv'
icPredictionsOutputFile = 'predictions/ExoCNN/ExoCNN.ic.extreme.90.unique.probabilities.csv'

evPredictions = fread(evPredictionsInputFile)
icPredictions = fread(icPredictionsInputFile)

RemoveDuplicatedAndSubstrings = function(predictions)
{
  ## Remove exactly duplicated seqs
  predictions = unique(predictions, by = 'seq')
  
  seqsArray = as.character(predictions$seq)
  names(seqsArray) = predictions$id
  
  seqsList = as.list(predictions$seq)
  
  ## Detect substrings
  isSubstring = mclapply(seqsList, function(aSeq)
    ifelse(sum(str_detect(seqsArray, aSeq), na.rm = T) > 1, TRUE, FALSE),
    mc.cores = detectCores() - 2)
  
  isSubstring = unlist(isSubstring)
  # table(isSubstring)
  
  ## Remove substrings
  predictions = predictions[!isSubstring]
  predictions  
}

evPredictions = RemoveDuplicatedAndSubstrings(evPredictions) # 503 seqs are removed
icPredictions = RemoveDuplicatedAndSubstrings(icPredictions) # 23 seqs are removed

fwrite(evPredictions, file = evPredictionsOutputFile)
fwrite(icPredictions, file = icPredictionsOutputFile)
