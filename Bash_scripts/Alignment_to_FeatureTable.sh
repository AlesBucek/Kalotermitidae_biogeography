#!/bin/bash

###set up input variables and paths to input data
#IMPORTANT: note on sample IDs: IDs of genomes should match between input alignment block, fasta headers of input genomes, and file names of input genome fastas
#INPORTANT: sample IDs MUST NOT include following strings and characters (without parentheses): "/", $feature variable as defined below, "END", "START"
#IMPORTANT: for each feature, the script should be run in separate working directory
#feature name (WARNING: must uniquely match the "product" name of target feature in feature table - either full match or partial match. E.g. "16S" is good match for "16S ribosomal RNA"))
feature=$1
#copy-pasted block of alignment from JalView (WARNING: must be named as follows: ${feature}.txt )
InBlock=$2
#feature teamplate which will be updated by feature coordinates (WARNING: prepare template for each feature)
FeatureTemplate=$3
#directory with original feature tables (WARNING: feature tables must have unique names matching the sequence names in fasta headers of ${InBlock} followed by sufix "_mitoscaf.fa.tbl" )
InTbls=$4
#directory containing files with fasta sequences of target genomes (WARNING: must be one genome per file, whole sequence on single line, file must have a name matching the sequence name in ${InBlock} followed by ".fasta" suffix)
GenomeDir=$5

#2) polish the copy-pasted block of multiple alignment:
#remove coordinates added by JalView to fasta header
sed 's/\/.*//' $InBlock > InBlock_polished 
#remove alignment gaps ("-")
sed -i '/^>/! s/-//g' InBlock_polished
# convert to single line-per-sequence format
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' InBlock_polished | tr "\t" "\n" > InBlock_polished_singleline
#split features into one sequence-per-file
awk '/^>/ {OUT=substr($0,2) ".fasta"}; OUT {print >OUT}' InBlock_polished_singleline

#3) find exact match for each feature sequence in target genome (WARNING: assumes single match, TODO: add check whether there is a single hit)
#loop through mt genomes (in "GenomeDir") to extract coordinates of features (in "FeatureDir"), input them into user-provided feature template, and output the filled-in feature templates into "${FeatureDir}/${FeatureTemplate}_${SampleID}_filled"
for GenomeFile in ${GenomeDir}/*.fasta
do SampleID=$(grep ">" $GenomeFile | sed 's/>//')
FeatureSeq=$(tail -n 1 ${SampleID}.fasta)
FeatureLen=$(echo $FeatureSeq | wc -c )
GenomeSeq=$(tail -n 1 $GenomeFile)
start=$(awk -v FeatureSeq="$FeatureSeq" 'i=index($0, FeatureSeq) {print i}' <(echo $GenomeSeq))
end=$(($start + $FeatureLen - 2)) 
#generate template only if feature length is >1 (i.e. feature is present)
if [ $FeatureLen -gt 1 ] 
then 
cp  ${FeatureTemplate} FeatureTemplate_${SampleID}_filled
sed -i "s/START/$start/" FeatureTemplate_${SampleID}_filled
sed -i "s/END/$end/" FeatureTemplate_${SampleID}_filled
fi
done

#4) remove existing feature annotations from feature tables
#directory to be built with feature tables from $InTbls after deleting the target feature
#remove 4 lines of feature table corresponding to the target feature
cp -r $InTbls FeatureTables_cleaned
for file in FeatureTables_cleaned/*.tbl ; do tac $file | sed -e "/${feature}/,+3d" | tac > ${file}2 ; done

#5) add features from filled-in feature template to cleaned feature tables
mkdir  FeatureTables_cleaned_updated
for file in FeatureTables_cleaned/*_mitoscaf.fa.tbl2 ; \
do ID=$(basename "$file" _mitoscaf.fa.tbl2) ; \
cat ${file} FeatureTemplate_${ID}_filled > FeatureTables_cleaned_updated/${ID}_mitoscaf.fa.tbl ; \
done
	#=> updated feature tables in FeatureTables_cleaned_updated/*_mitoscaf.fa.tbl
