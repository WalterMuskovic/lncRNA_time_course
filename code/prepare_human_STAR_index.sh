#!/bin/bash

## Prepare GRCh38 STAR index for 2-pass alignment.

# Create first STAR index. This will be rebuilt after first-pass mapping of all totalRNA samples.
mkdir $data_dir/STAR_index
STAR --runMode genomeGenerate --genomeFastaFiles $data_dir/genome/GRCh38_spiked.fa --sjdbGTFfile $data_dir/annotation/GRCh38_spiked.gtf --sjdbOverhang 124 --runThreadN 10 --genomeDir $data_dir/STAR_index --outFileNamePrefix STAR_GRCh38

# Align all samples, combine the *SJ.out.tab and rebuild the STAR genome index with novel junctions included
mkdir $data_dir/temp_STAR
cd $data_dir/temp_STAR
for time_point in {000..400..10}
do
STAR --genomeDir $data_dir/STAR_index/ --readFilesIn "${data_dir}/trimmed_fastq/trimmed_T98G_${time_point}_min_R1.FASTQ.gz" ".data/trimmed_fastq/trimmed_T98G_${time_point}_min_R2.FASTQ.gz" \
    --readFilesCommand zcat --runThreadN 12 --genomeLoad NoSharedMemory \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile $data_dir/temp_STAR/COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
    --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
    --limitBAMsortRAM 60000000000 --outFileNamePrefix "${data_dir}/temp_STAR/T98G_${time_point}_min" \
    --sjdbGTFfile $data_dir/annotation/GRCh38_spiked.gtf \
    --alignEndsType Local
done

# rebuild STAR index
# concatenate the SJ.out.tab files produced by all alignments
cat $data_dir/temp_STAR/*SJ.out.tab > $data_dir/STAR_index/combined.SJ.out.tab 

# Re-build STAR genome index, supplying combined junctions file combined.SJ.out.tab from first pass of all time points
STAR --runMode genomeGenerate --genomeFastaFiles $data_dir/genome/GRCh38_spiked.fa --sjdbGTFfile $data_dir/annotation/GRCh38_spiked.gtf --sjdbFileChrStartEnd $data_dir/STAR_index/combined.SJ.out.tab --limitSjdbInsertNsj 1500000 --sjdbOverhang 124 --runThreadN 10 --genomeDir $data_dir/STAR_index --outFileNamePrefix STAR_T98G
# Clear temporary directory
rm -rf $data_dir/temp_STAR

echo "Finished building STAR GRCh38 genome index for two-pass alignment of T98G samples"
