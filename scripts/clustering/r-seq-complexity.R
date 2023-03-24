setwd('/media/pgdrive/sharif/exosomians/')

options(stringsAsFactors = F)
options(scipen = 999)

library(data.table)

evPredictionsFile = 'predictions/ExoCNN/ExoCNN.ev.extreme.90.unique.probabilities.csv'

evx = fread(evPredictionsFile)

