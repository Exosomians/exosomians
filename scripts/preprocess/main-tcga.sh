#!/bin/bash
set -e

sampleType=$1
dataDir=$2

export sampleType=$sampleType

# @ TODO
### THRESHOLD TO BE CONSIDERED:
# THRESHOLD1: MINIMUM VALID READS COVERAGE
# THRESHOLD2: MINIMUM VALID RNA LENGTH
# THRESHOLD3: MAXIMUM VALID RNA LENGTH

THREADS=`nproc`
THREADS=`expr $THREADS - 2`

minLength=18
maxLength=500

# for f in $dataDir/*.bam; do
#  samtools sort -@ $THREADS -o ${f/.bam/.sorted.bam} $f;
# done

# samtools merge -@ $THREADS ${sampleType}.all.sorted.bam $dataDir/*.sorted.bam
 
### THRESHOLD1
## check number of all reads to find a reasonable coverage threshold proportional to the total read numbers
## final threshold is max of 5 and the proportional threshold (2/3 of milion reads of input sample's lib size)
## samtools flagstat all_combined.sort.bam

# samtools view ${sampleType}.all.sorted.bam | wc -l > library.size

# python3 <<CODE
# from math import ceil
# with open('library.size', 'r') as f:
# 	libSize = int(f.read())
# minCoverage = ceil(libSize / 1e6 * (2/3))
# minCoverage = max(minCoverage, 5)
# with open('minimum.coverage', 'w') as f:
# 	f.write(str(minCoverage))
# CODE

# minCoverage=`cat minimum.coverage`
# echo Minimum Valid Coverage is $minCoverage
minCoverage=20

## extract regions with minimum coverage threshold for each negative and positive strands separately
cat ${sampleType}.pos_strand_regions.bed | \
	awk -v minCov="$minCoverage" '$4>minCov' > ${sampleType}.pos_strand_covered_regions.bed

cat ${sampleType}.neg_strand_regions.bed | \
	awk -v minCov="$minCoverage" '$4>minCov' > ${sampleType}.neg_strand_covered_regions.bed

cat ic.pos_strand_covered_regions.bed >> ${sampleType}.pos_strand_covered_regions.bed
cat ic.neg_strand_covered_regions.bed >> ${sampleType}.neg_strand_covered_regions.bed

bedtools sort -i ${sampleType}.pos_strand_covered_regions.bed > ${sampleType}.pos_strand_covered_regions.sorted.bed
bedtools sort -i ${sampleType}.neg_strand_covered_regions.bed > ${sampleType}.neg_strand_covered_regions.sorted.bed

rm ${sampleType}.pos_strand_covered_regions.bed 
rm ${sampleType}.neg_strand_covered_regions.bed

mv ${sampleType}.pos_strand_covered_regions.sorted.bed ${sampleType}.pos_strand_covered_regions.bed 
mv ${sampleType}.neg_strand_covered_regions.sorted.bed ${sampleType}.neg_strand_covered_regions.bed 

# filter regions by max 500 and min 18 length
bedtools merge -i ${sampleType}.pos_strand_covered_regions.bed | \
	awk -v maxLen="$maxLength" '($3-$2)<=maxLen' | \
	awk -v minLen="$minLength" '($3-$2)>=minLen' > ${sampleType}.pos_inrange_covered_regions.bed

bedtools merge -i ${sampleType}.neg_strand_covered_regions.bed | \
	awk -v maxLen="$maxLength" '($3-$2)<=maxLen' | \
	awk -v minLen="$minLength" '($3-$2)>=minLen' > ${sampleType}.neg_inrange_covered_regions.bed

# python3 <<CODE
# with open('pos_inrange_covered_regions.bed', 'r') as istrp:
#   with open('pos_inrange_covered_regions_stranded.bed', 'w') as ostrp:
#     for line in istrp:
#       line = line.rstrip('\n') + '\t.\t0\t+'
#       print(line, file=ostrp)

# with open('neg_inrange_covered_regions.bed', 'r') as istrn:
#   with open('neg_inrange_covered_regions_stranded.bed', 'w') as ostrn:
#     for line in istrn:
#       line = line.rstrip('\n') + '\t.\t0\t-'
#       print(line, file=ostrn)
# CODE


python3 <<CODE
import os
sampleType = str(os.environ["sampleType"])

with open(f'{sampleType}.pos_inrange_covered_regions.bed', 'r') as istrp:
  with open(f'{sampleType}.pos_inrange_covered_regions_stranded.bed', 'w') as ostrp:
    for line in istrp:
      line = line.rstrip('\n') + '\t.\t0\t+'
      print(line, file=ostrp)

with open(f'{sampleType}.neg_inrange_covered_regions.bed', 'r') as istrn:
  with open(f'{sampleType}.neg_inrange_covered_regions_stranded.bed', 'w') as ostrn:
    for line in istrn:
      line = line.rstrip('\n') + '\t.\t0\t-'
      print(line, file=ostrn)
CODE



cat ${sampleType}.pos_inrange_covered_regions_stranded.bed ${sampleType}.neg_inrange_covered_regions_stranded.bed | \
	bedtools sort -i stdin > ${sampleType}.final.all_inrange_covered_regions_stranded.bed
	# sort -k 1,1 2,2n > all_inrange_covered_regions_stranded.bed

cp ${sampleType}.final.all_inrange_covered_regions_stranded.bed ic.plus.tcga.final.all_inrange_covered_regions_stranded.bed

#rm \
#	library.size \
#	minimum.coverage \
#	${sampleType}.pos_inrange_covered_regions.bed \
#	${sampleType}.pos_strand_covered_regions.bed \
#	${sampleType}.neg_inrange_covered_regions.bed \
#	${sampleType}.neg_strand_covered_regions.bed \
#	${sampleType}.pos_inrange_covered_regions_stranded.bed \
#	${sampleType}.neg_inrange_covered_regions_stranded.bed

