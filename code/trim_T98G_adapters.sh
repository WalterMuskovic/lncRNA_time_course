#!/bin/bash

mkdir $data_dir/trimmed_fastq
cd $data_dir/fastw

# Using cutadapt 1.11 with Python 2.7.5
for time_point in {000..400..10}
do 
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  --minimum-length 15 \
  -o "${data_dir}/trimmed_fastq/trimmed_T98G_${time_point}_min_R1.FASTQ.gz" \
  -p "${data_dir}/trimmed_fastq/trimmed_T98G_${time_point}_min_R2.FASTQ.gz" \
  "T98G_${time_point}_min_R1.FASTQ.gz" "T98G_${time_point}_min_R2.FASTQ.gz"
done

# Create md5 checksums for trimmed fastq files
md5sum $data_dir/trimmed_fastq/*.fastq.gz > $data_dir/trimmed_fastq/trimmed_fastq_md5.chk
