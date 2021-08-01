#!/bin/bash

# Define data directory
data_dir=$1
data_track_dir=$2

cd $data_dir/fastq
cp $data_track_dir/T98G_md5_checklist.chk $data_dir/fastq/T98G_md5_checklist.chk
md5sum -c T98G_md5_checklist.chk

# Create file with number of lines in each fastq file (divide by 4 to get number of reads per sample)
touch T98G_counts.txt
for time_point in {000..400..10}
do
zcat "T98G_${time_point}_min_R1.fastq.gz"  | wc -l >> T98G_counts.txt
done

cp $data_dir/fastq/T98G_counts.txt $data_track_dir/T98G_counts.txt

echo "Finished checking md5sums for T98G fastq files"
