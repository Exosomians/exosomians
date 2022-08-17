#!/bin/bash
set -eo pipefail

memeFiles=$1
while read memeFile;
do
	sed -n '/BL /,/\/\//p' $memeFile > ${memeFile/.txt/_motifs.txt};
	csplit -b "%04d.txt" ${memeFile/.txt/_motifs.txt} "/BL/" "{*}";
	sed -i '1d;$d' xx*.txt && sed -i 's/ (.*//g' xx*.txt;
	rm xx0000.txt
	for ids in xx**.txt; 
	do 
		grep -f $ids /media/pgdrive/sharif/exosomians/predictions/ExoCNN/ExoCNN.ev.extreme.90.probabilities.csv > ${ids/.txt/.seqs.ids};
		sed -i '1 i\,id,seq,no,yes' ${ids/.txt/.seqs.ids}; 
	done; 
	
	dirName=`dirname $memeFile`;
	dirName=${dirName/\//-};
	for xx in xx**; 
	do 
		
		mv $xx ${xx/xx/${dirName}_};
		# rm xx*.txt && rm xx0000.*;
	done;
	rm *.txt
done < $memeFiles
