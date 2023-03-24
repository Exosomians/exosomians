#!/bin/bash
set -e

homerDir=$1

# findMotifs.pl targets.fa fasta motifResults/ -fasta background.fa 

# $homerDir/cpp/homer2 denovo \
# 	-i ev.extreme.dna.fasta \
# 	-b ic.extreme.dna.fasta \
#	-len 10 \
#	-p 32 \
#	-cache 16000 \
#	> $homerDir/homerData/homer2.denovo.10.results

findMotifs.pl \
	ev.extreme.dna.fasta \
	fasta \
	$homerDir/homerData/final/ \
	-fasta ic.extreme.dna.fasta \
	-p 32 \
	-cache 16000 \
	-len 4,6,8 \
	-S 50 \
	-rna

