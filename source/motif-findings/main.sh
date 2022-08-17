#!/bin/bash

# @ TODO
### THRESHOLD TO BE CONSIDERED:
# THRESHOLD1: MINIMUM VALID READS COVERAGE
# THRESHOLD2: MINIMUM VALID RNA LENGTH
# THRESHOLD3: MAXIMUM VALID RNA LENGTH

## Sort all IC and EV bams (Hani's data are already sorted)
# for f in IC/*.bam; do
#   out="$f.sort.bam"
#  samtools sort -o $out $f
# done
# for f in EV/*.bam; do
#   out="$f.sort.bam"
#   samtools sort -o $out $f
# done

#### MERGING ALL IC SAMPLES ####

icDir=$1

samtools merge all_ICs_combined.sort.bam $icDir/*.bam

### THRESHOLD1
## check number of all reads to find a reasonable coverage threshold proportional to the total read numbers
## final threshold is max of 5 and the proportional threshold (2/3 of milion reads of input sample's lib size)
# samtools flagstat all_ICs_combined.sort.bam

# extract regions with minimum coverage threshold for each negative and positive strands separately
bedtools genomecov -ibam all_ICs_combined.sort.bam -bg -strand + | \
	awk '$4>20' > pos_strand_covered_regions.bed
bedtools genomecov -ibam all_ICs_combined.sort.bam -bg -strand - | \
	awk '$4>20' > neg_strand_covered_regions.bed

# filter regions by max 500 and min 18 length
bedtools merge -i pos_strand_covered_regions.bed | \
	awk '($3-$2)<=500' | \
	awk '($3-$2)>=18' > pos_inrange_covered_regions.bed

bedtools merge -i neg_strand_covered_regions.bed | \
	awk '($3-$2)<=500' | \
	awk '($3-$2)>=18' > neg_inrange_covered_regions.bed

python3 <<CODE
	with open(pos_inrange_covered_regions.bed, 'r') as istrp:
	  with open(pos_inrange_covered_regions_stranded.bed, 'w') as ostrp:
	    for line in istrp:
	      line = line.rstrip('\n') + '\t.\t0\t+'
	      print(line, file=ostrp)

	with open(neg_inrange_covered_regions.bed, 'r') as istrn:
	  with open(neg_inrange_covered_regions_stranded.bed, 'w') as ostrn:
	    for line in istrn:
	      line = line.rstrip('\n') + '\t.\t0\t-'
	      print(line, file=ostrn)
CODE


cat pos_inrange_covered_regions_stranded.bed neg_inrange_covered_regions_stranded.bed | \
	bedtools sort -i stdin > all_inrange_covered_regions_stranded.bed
	# sort -k 1,1 2,2n > all_inrange_covered_regions_stranded.bed