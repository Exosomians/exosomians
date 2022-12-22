#! /usr/bin/env python

#shuffleSequenceUsingAltschulErikson.py
#P. Clote, Oct 2003

#------------------------------------------------------------------
#Input RNAs in FASTA file, and compute NUM many shufflings of RNA sequence
#using Altman-Erikson randomly shuffled dinucleotide method.
#------------------------------------------------------------------

PRINT   = 0
LINELEN = 70

import sys, os, string
import numpy as np

from exopy.generator._erikson_dinucl_shuffle import dinuclShuffle


def generate_random_sequence(sequences, NUM):
  random_sequences = []
  for sequence in sequences.tolist():
    for i in range(NUM):
      shuffledSeq = dinuclShuffle(sequence) 
      random_sequences.append(shuffledSeq)
  return np.array(random_sequences)

  
# if __name__ == '__main__':  
#   if len(sys.argv) < 3:
#      print("Usage: %s RNAs.faa NUM" %  sys.argv[0])
#      text = """
#             1) RNA sequence
#             2) NUM is number of random sequences to generate by
#                shuffling the dinucleotides of RNAs input
#      Script to compute Altman-Erikson randomly shuffled dinucleotides.
#             """
#      print(text)
#      sys.exit(1)
#   sequence = sys.argv[1]
#   NUM      = int(sys.argv[2])
#   generate(sequence, NUM)
  