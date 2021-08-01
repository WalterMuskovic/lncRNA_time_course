#!/bin/bash

# Genome sequence
mkdir $data_dir/genome
wget -P $data_dir/genome ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
# Checksum
cd $data_dir/genome
md5sum GRCh38.primary_assembly.genome.fa.gz | grep e319d85f95bb780bcd255d97ab20b118
# Unzip
gunzip GRCh38.primary_assembly.genome.fa.gz

# Sequin sequences
wget https://s3.amazonaws.com/sequins/annotations/chrIS.v2.fa
# Checksum
md5sum chrIS.v2.fa | grep 9bea7e851721f6360d09bceccb8c4371

# Combine sequences
cat GRCh38.primary_assembly.genome.fa chrIS.v2.fa > GRCh38_spiked.fa
# Verify chromosomes
grep chr GRCh38_spiked.fa
# checks out

# Remove individual sequence files
rm GRCh38.primary_assembly.genome.fa chrIS.v2.fa

echo "Finished downloading human genome sequence and sequin sequences"
