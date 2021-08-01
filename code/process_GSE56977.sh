#!/bin/bash

# In addition to the human glioblastoma T98G time course data we are going to access a publicly available
# mouse time course dataset described in [Rabani et al, Cell,  2014](https://doi.org/10.1016/j.cell.2014.11.015).
# The data is available using this accession [GSE56977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56977).
# The dataset describes mouse dendritic cells responding to LPS at 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165
# and 180 min intervals post-treatment, with one total RNA-seq sample per time point.

# Define directories
data_dir=$1
data_track_dir=$2

mkdir $data_dir/mouse_temp
cd $data_dir/mouse_temp

##The following code will retrieve and process the GSE56977 time series data (SRA > FASTQ > (QC + trimming) > BAM).
# Prefetch Rabani et al SRA files
for accession in {1258347..1258359}
do
  prefetch -O "SRR${accession}" &
done
wait

# Extract fastq files
for accession in {1258347..1258359}
do
  fastq-dump.2 --gzip --split-3 -O ./ "SRR${accession}.sra" &
done
wait

# Check how many lines are present in sample_name_R1.fastq.gz and sample_name_R2.fastq.gz files - should be equal.
# Also check how many reads are in the sample_name.fastq.gz file - if present.
for accession in {1258347..1258359}
do
  zcat "SRR${accession}_1.FASTQ.gz" | wc -l
  zcat "SRR${accession}_2.FASTQ.gz" | wc -l
  zcat "SRR${accession}.FASTQ.gz" | wc -l
done

# Create file with number of lines in each fastq file (divide by 4 to get number of reads per sample)
touch $data_track_dir/mouse_LPS_counts.txt
for accession in {1258347..1258359}
do
zcat "SRR${accession}_1.fastq.gz"  | wc -l >> $data_track_dir/mouse_LPS_counts.txt
done

# Run fastQC on samples
mkdir fastqc
for accession in {1258347..1258359}
do
  fastqc --threads 10 --outdir fastqc/Rabani "SRR${accession}_1.FASTQ.gz" "SRR${accession}_2.FASTQ.gz" &
done
wait

# Trim Illumina TruSeq adapter sequences Using cutadapt 1.11 with Python 2.7.5
for accession in {1258347..1258359}
do 
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  --minimum-length 15 \
  -o "trimmed_SRR${accession}_1.FASTQ.gz" \
  -p "trimmed_SRR${accession}_2.FASTQ.gz" \
  "SRR${accession}_1.FASTQ.gz" "SRR${accession}_2.FASTQ.gz"
done

# Give files meaningful names and move to the trimmed_fastq directory
mv trimmed_SRR1258347_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_000_min_R1.fastq.gz
mv trimmed_SRR1258348_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_015_min_R1.fastq.gz
mv trimmed_SRR1258349_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_030_min_R1.fastq.gz
mv trimmed_SRR1258350_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_045_min_R1.fastq.gz
mv trimmed_SRR1258351_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_060_min_R1.fastq.gz
mv trimmed_SRR1258352_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_075_min_R1.fastq.gz
mv trimmed_SRR1258353_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_090_min_R1.fastq.gz
mv trimmed_SRR1258354_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_105_min_R1.fastq.gz
mv trimmed_SRR1258355_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_120_min_R1.fastq.gz
mv trimmed_SRR1258356_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_135_min_R1.fastq.gz
mv trimmed_SRR1258357_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_150_min_R1.fastq.gz
mv trimmed_SRR1258358_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_165_min_R1.fastq.gz
mv trimmed_SRR1258359_1.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_180_min_R1.fastq.gz

mv trimmed_SRR1258347_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_000_min_R2.fastq.gz
mv trimmed_SRR1258348_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_015_min_R2.fastq.gz
mv trimmed_SRR1258349_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_030_min_R2.fastq.gz
mv trimmed_SRR1258350_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_045_min_R2.fastq.gz
mv trimmed_SRR1258351_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_060_min_R2.fastq.gz
mv trimmed_SRR1258352_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_075_min_R2.fastq.gz
mv trimmed_SRR1258353_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_090_min_R2.fastq.gz
mv trimmed_SRR1258354_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_105_min_R2.fastq.gz
mv trimmed_SRR1258355_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_120_min_R2.fastq.gz
mv trimmed_SRR1258356_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_135_min_R2.fastq.gz
mv trimmed_SRR1258357_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_150_min_R2.fastq.gz
mv trimmed_SRR1258358_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_165_min_R2.fastq.gz
mv trimmed_SRR1258359_2.fastq.gz ${data_dir}/trimmed_fastq/Rabani_total_180_min_R2.fastq.gz

# Create STAR genome index for GRCm38
mkdir $data_dir/STAR_index_mouse_101bp
STAR --runMode genomeGenerate --genomeFastaFiles $data_dir/genome_mouse/GRCm38.primary_assembly.genome.fa --sjdbGTFfile $data_dir/annotation_mouse/gencode.vM21.annotation.gtf --sjdbOverhang 100 --runThreadN 10 --genomeDir $data_dir/STAR_index_mouse_101bp --outFileNamePrefix STAR_GRCm38_101bp

# Align samples using the newly created genome index
for time_point in total_000 total_015 total_030 total_045 total_060 total_075 total_090 total_105 total_120 total_135 total_150 total_165 total_180
do
# Align
cd $data_dir/alignments
STAR --genomeDir $data_dir/STAR_index_mouse_101bp --readFilesIn "${data_dir}/trimmed_fastq/Rabani_${time_point}_min_R1.FASTQ.gz" "${data_dir}/trimmed_fastq/Rabani_${time_point}_min_R2.FASTQ.gz" \
    --readFilesCommand zcat --runThreadN 10 --genomeLoad NoSharedMemory \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
    --quantMode GeneCounts --sjdbScore 1 \
    --limitBAMsortRAM 60000000000 --outFileNamePrefix "${data_dir}/alignments/Rabani_${time_point}_min" \
    --sjdbGTFfile $data_dir/annotation_mouse/gencode.vM21.annotation.gtf \
    --alignEndsType Local

  # Create BAM index file
  samtools index "Rabani_${time_point}_minAligned.sortedByCoord.out.bam"
done

# clean up files that are no longer needed
rm -rf $data_dir/mouse_temp
