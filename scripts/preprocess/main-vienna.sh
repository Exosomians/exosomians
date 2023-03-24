#!/bin/bash
# wget https://www.tbi.univie.ac.at/RNA/download/ubuntu/ubuntu_18_04/viennarna_2.4.14-1_amd64.deb
# sudo apt install libgsl23 libgslcblas0
# sudo dpkg -i viennarna_2.4.14-1_amd64.deb
# sudo apt install -f

JOBS=`nproc`
JOBS=`expr $JOBS - 2`
SEQUENCES_FASTA=$1
# DOTBRACKET_FASTA=$2

RNAfold  --noPS --jobs=$JOBS -i $SEQUENCES_FASTA | \
       grep -e ^\> -e ^\( -e ^\[.] \
       > ${SEQUENCES_FASTA/.fasta/.dotbracket.fasta}
