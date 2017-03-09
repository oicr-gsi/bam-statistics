#Run Samtools depth on your BAM files and save output as depth.txt 

/PATH/TO/SAMTOOLS/samtools depth -b PATH/TO/BED_FILE -f /PATH/TO/FILE/CONTAINING/LIST/OF/BAM_FILES > depth.txt

# Change locations/position in each interval to index (0 to 119) using python script changeLocationToIndex.py 
cut -f2 depth.txt > positions.txt
python /PATH/TO/SCRIPT/changeLocationToIndex.py positions.txt > index.txt

# Get coverage values for each BAM
# I had 6 different BAM files so my depth output contains 6 different columns corresponding to each BAM
cut -f3,4,5,6,7,8 depth.txt > cov.txt

# [Optional] Use python script normalize.py to normalize coverage of all position in each BAM
python /PATH/TO/SCRIPT/normalize.py cov.txt > norm_cov.txt

# Concatenate index file and normalized coverage 
paste  index.txt norm_cov.txt | column -s $'\t' -t > f7.txt 

# Calculate rows with similar indices and average of coverage at each index for each BAM using bash script avgAllColumns.py
bash PATH/TO/SCRIPT/avgAllColumns.sh > FILENAME.csv

# Change file format from tsv to csv
sed -i 's/\t/,/g' FILENAME.csv 

# [optional] Add Approriate Header to CSV file
echo 'distance,Dataset_1,Dataset_2,Dataset_3,Dataset_4,Dataset_5,Dataset_6' > header.csv 
cat header.csv FILENAME.csv > FINAL_OUTPUT_FILE.csv

# Clean up (remove all temp files)
rm header.csv FILENAME.csv f7.txt norm_cov.txt cov.txt index.txt position.txt depth.txt
