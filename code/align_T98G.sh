#!/bin/bash

mkdir $data_dir/alignments
cd $data_dir/alignments

for time_point in {000..400..10}
do
# Carry out STAR 2nd-pass paired-end read mapping
STAR --genomeDir $data_dir/STAR_index/ --readFilesIn "${data_dir}/trimmed_fastq/trimmed_T98G_${time_point}_min_R1.FASTQ.gz"  "${data_dir}/trimmed_fastq/trimmed_T98G_${time_point}_min_R2.FASTQ.gz"  \
    --readFilesCommand zcat --runThreadN 12 --genomeLoad NoSharedMemory \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile $data_dir/alignments/COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
    --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 --limitSjdbInsertNsj 1500000 \
    --limitBAMsortRAM 60000000000 --outFileNamePrefix "${data_dir}/alignments/T98G_${time_point}_min" \
    --sjdbGTFfile $data_dir/annotation/GRCh38_spiked.gtf \
    --alignEndsType Local

# Create BAM index file with samtools-1.3.1
samtools index "T98G_${time_point}_minAligned.sortedByCoord.out.bam"
done

# Create md5 checksums for genome- and transcriptome-aligned bams
md5sum *_minAligned.sortedByCoord.out.bam > genome_aligned_bams_md5.chk
md5sum *_minAligned.toTranscriptome.out.bam > transcriptome_aligned_bams_md5.chk

echo "Finished aligning T98G time series samples to GRCh38"
