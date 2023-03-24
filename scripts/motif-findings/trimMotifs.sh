#!/bin/bash
set -eo pipefail

for ids in *.seqs.ids;
do
	numOfMotifs=`wc -l $ids | cut -f1 -d' '`;
	# echo $numOfMotifs
	if [ $numOfMotifs -lt 31 ]
	then
		rm $ids
		echo "$ids removed!"
	fi
done;
