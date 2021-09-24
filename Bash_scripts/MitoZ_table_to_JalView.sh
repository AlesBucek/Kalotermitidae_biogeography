#!/bin/bash

#variables and dirs
tblInDir=$1
gffOutDir=$2
mkdir ${gffOutDir}

#run loop
for file in ${tblInDir}/*_mitoscaf.fa.tbl
do while read line
do ID=$(basename "$file" _mitoscaf.fa.tbl)
feature=$(echo $line | cut -f3 -d" ")
#strip ">" and "<" from  coordinates of incomplete features
start=$(echo $line | cut -f1 -d" " | sed 's/>//' | sed 's/<//')
stop=$(echo $line | cut -f2 -d" " | sed 's/>//' | sed 's/<//')
#add frame (+ or -) by checking if start coordinate is smaller or larger than stop coordinate
#order coordinates so lower number is before higher number (gff format requirement)
let "res = $start - $stop"
if [ $res -lt 0 ]
then coordLow=$start
coordHigh=$stop
frame="+"
else coordLow=$stop
coordHigh=$start
frame="-"
fi
awk -v ID="$ID" -v frame="$frame" -v feature="$feature" -v coordLow="$coordLow" -v coordHigh="$coordHigh" -v OFS='\t' 'BEGIN{print ID, feature, feature, coordLow, coordHigh, ".", frame, "."}' >> ${gffOutDir}/${ID}.gff
#input into loop pre-processed feature tables (only lines starting with number or ">number" or "<number" ~ feature lines including incomplete features)
done < <(grep "^>[0-9]\|^<[0-9]\|^[0-9]" $file | grep -v "gene"  )
done

#concatenate GFFs and add header
cat <(echo "GFF") ${gffOutDir}/*.gff > ${gffOutDir}.txt
