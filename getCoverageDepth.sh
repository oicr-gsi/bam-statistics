#!/bin/bash

## This script takes in 4 different files.  The first one is a text file containing a list of BAM files,
## one per line.  The second file is the BED file of off target regions on the upstream boundaries of some 
## target intervals. The third file is the BED file of off target regions on the downstream boundaries of 
## the target intervals.  Finally the last file is the genome reference file.

## Usage: 
## bash getCoverageDepth.sh /PATH/TO/BAM_List_File \
##      PATH/TO/BED_File_Upstream \
##      PATH/TO/BED_File_Downstream \
##      PATH/TO/Reference_File


# ADD YOUR PATH TO SAMTOOLS HERE!!!!!
SAMTOOLS=/u/sbundhoo/tools/samtools/

BAM_LIST_FILE="$1"
BED_FILE_UPSTREAM="$2"
BED_FILE_DOWNSTREAM="$3"
REFERENCE_FILE="$4"

if [ ! -f $BAM_LIST_FILE ]; then 
  echo "File $BAM_LIST_FILE does not exist" 
  exit
fi

if [ ! -f $BED_FILE_UPSTREAM ]; then
  echo "File $BED_FILE_UPSTREAM does not exist" 
  exit
fi

if [ ! -f $BED_FILE_DOWNSTREAM ]; then
  echo "File $BED_FILE_DOWNSTREAM does not exist"
  exit
fi

if [ ! -f $REFERENCE_FILE ]; then
  echo "File $REFERENCE_FILE does not exist"
  exit
fi


## Calculate coverage depth on upstream target boundary regions
$SAMTOOLS/samtools depth -a -b $BED_FILE_UPSTREAM -f $BAM_LIST_FILE --reference $REFERENCE_FILE > ./upstream_temp1.txt

cut -f2 upstream_temp1.txt | python ./changePositionToIndexUpstream.py > upstream_temp2.txt

## Note that my BAM_LIST_FILE contains 6 different BAM files, that's we have 6 columns(3,4,5,6,7,8) 
## corresponding to the coverage depth of each BAM:

cut -f3,4,5,6,7,8 upstream_temp1.txt > upstream_temp3.txt 

## coverage depths can be normalized using normalize.py 

paste upstream_temp2.txt upstream_temp3.txt | column -s $'\t' -t > upstream_temp4.txt 



## Calculate coverage depth on downstream target boundary regions
$SAMTOOLS/samtools depth -a -b $BED_FILE_DOWNSTREAM -f $BAM_LIST_FILE --reference $REFERENCE_FILE > ./downstream_temp1.txt


cut -f2 downstream_temp1.txt | python ./changePositionToIndexDownstream.py > downstream_temp2.txt

## Note that my BAM_LIST_FILE contains 6 different BAM files, that's we have 6 columns(3,4,5,6,7,8) 
## corresponding to the coverage depth of each BAM:

cut -f3,4,5,6,7,8 downstream_temp1.txt > downstream_temp3.txt 

## coverage depths can be normalized using normalize.py script 

paste downstream_temp2.txt downstream_temp3.txt | column -s $'\t' -t > downstream_temp4.txt 

cat upstream_temp4.txt downstream_temp4.txt > combined_cov.txt

## calculate average of coverage depths at same indices

awk '
  {
  c[$1]++; 
  for (i=2;i<=NF;i++) {
    s[$1"."i]+=$i};
  } 
  END {
    for (k in c) {
      printf "%s\t", k; 
      for(i=2;i<NF;i++) printf "%.3f\t", s[k"."i]/c[k]; 
      printf "%.3f\n", s[k"."NF]/c[k];
#      printf "%s\n", c[k];
    }
  }' combined_cov.txt > cov.txt 


# Change file format from tsv to csv
sed -i 's/\t/,/g' cov.txt 

# [optional] Add Approriate Header to CSV file
echo 'distance,Dataset_1,Dataset_2,Dataset_3,Dataset_4,Dataset_5,Dataset_6' > header.csv 
cat header.csv cov.txt > cov.csv

##Clean up
rm upstream_temp* downstream_temp* combined_cov.txt header.csv cov.txt 

## Coverage by distance can be plotted using R
