#!/bin/bash

cd $data_dir/alignments

# Merge BAM files
samtools merge --threads 64 t98g_merged.bam -b $data_dir/t98g_bam_file_list.txt 
samtools index t98g_merged.bam
samtools merge --threads 64 Rabani_merged.bam -b $data_dir/Rabani_bam_file_list.txt 
samtools index Rabani_merged.bam

# Run StringTie
stringtie t98g_merged.bam -G $data_dir/annotation/GRCh38_spiked.gtf -o $data_dir/annotation/T98G_stringtie.gtf -j 1 -a 10 -A $data_dir/annotation/T98G_gene_abundance.tab -m 200 -l T98G_STRG -p 64 -v --rf -x 'chrM,chrIS'
stringtie Rabani_merged.bam -G $data_dir/annotation/gencode.vM21.annotation.gtf -o $data_dir/annotation/Rabani_stringtie.gtf -j 1 -a 10 -A $data_dir/annotation/Rabani_gene_abundance.tab -m 200 -l Rabani_STRG -p 64 -v --rf -x 'chrM'

# Remove small number of features that are unstranded (strand=".") - causes problems later on when we try to build a TxDb
awk '$7 != "."' $data_dir/annotation/T98G_stringtie.gtf > $data_dir/annotation/T98G_stringtie_stranded.gtf
awk '$7 != "."' $data_dir/annotation/Rabani_stringtie.gtf > $data_dir/annotation/Rabani_stringtie_stranded.gtf
