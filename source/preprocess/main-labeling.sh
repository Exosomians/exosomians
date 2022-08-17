#!/bin/bash
set -e

icBed=$1
evBed=$2

minSubstractResidueLength=18

bedtools subtract \
	-a $icBed \
	-b $evBed \
	-s | \
	awk -v minSubResLen=$minSubstractResidueLength '($3-$2)>=minSubResLen' \
	> ic.final.bed

mv ${evBed} ev.final.bed