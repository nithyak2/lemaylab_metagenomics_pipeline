#!/bin/bash
# prompt entry of fastq file path
read -p "Enter path to fastq files: " path
# prompt entry of output paths
read -p "Enter output filename: " output

# compute the md5 checksums for each file
for file in $path*.gz; do
	md5sum $file >> $output
	sed -i "s|$path||g" $output
done

# remove the file path string in the sequence name header

exit 0	
