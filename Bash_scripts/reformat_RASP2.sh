#!/bin/bash

#define input and output dirs and files
indir=$1
outdir=$2
inlist=$3
mkdir $outdir

#run reformating in loop
while read line ;\
do infile=$indir/$line ;\
fbasename=$(basename "$line" .txt) ;\
while read line ;\
do node=$(echo $line | cut -d" " -f2) ;\
tip1num=$(echo $line | grep -o -P '(?<=terminals ).*(?=-)') ;\
tip1ID=$(awk '/\[TAXON\]/,/\[TREE\]/' $infile | grep -v "\[TAXON\]" | grep -v "\[TREE\]" | grep "^$tip1num[[:space:]]" | grep -o -P '(?<=\t).*(?=\t)');\
tip2num=$(echo $line | grep -o -P '(?<=-).*(?=\))') ;\
tip2ID=$(awk '/\[TAXON\]/,/\[TREE\]/' $infile | grep -v "\[TAXON\]" | grep -v "\[TREE\]" | grep "^$tip2num[[:space:]]" | grep -o -P '(?<=\t).*(?=\t)');\
probsA=$(echo $line | grep -o -P '(?<=A )[0-9\.]*') ;\
probsB=$(echo $line | grep -o -P '(?<=B )[0-9\.]*') ;\
probsC=$(echo $line | grep -o -P '(?<=C )[0-9\.]*') ;\
probsD=$(echo $line | grep -o -P '(?<=D )[0-9\.]*') ;\
probsE=$(echo $line | grep -o -P '(?<=E )[0-9\.]*') ;\
probsF=$(echo $line | grep -o -P '(?<=F )[0-9\.]*') ;\
probsG=$(echo $line | grep -o -P '(?<=G )[0-9\.]*') ;\
probsH=$(echo $line | grep -o -P '(?<=H )[0-9\.]*') ;\
probsI=$(echo $line | grep -o -P '(?<=I )[0-9\.]*') ;\
probsJ=$(echo $line | grep -o -P '(?<=J )[0-9\.]*') ;\
echo -e $node '\t' $tip1ID '\t' $tip2ID '\t' $probsA '\t' $probsB '\t' $probsC '\t' $probsD '\t' $probsE '\t' $probsF '\t' $probsG '\t' $probsH '\t' $probsI '\t' $probsJ  >> $outdir/${fbasename}.txt ;\
done < <(awk '/Result of combined:/,/Result of run 1:/' $infile | grep -v "Result of combined:" | grep -v "Result of run 1:") ;\
sed -i '1s;^;node\ttip1ID\ttip2ID\tAfrotropical\tOriental\tAustralian\tNeotropical\tOceanian\tMadagascan\tSaharo_Arabian\tNearctic\tSino_Japanese\tPalearctic\n;' $outdir/${fbasename}.txt ;\
done < "$inlist"

