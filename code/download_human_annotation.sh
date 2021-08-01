#!/bin/bash

# Genome annotation
mkdir $data_dir/annotation
cd $data_dir/annotation
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
# Checksum
md5sum gencode.v29.annotation.gtf.gz | grep 949dca487934459e6338090040e77628  
# Unzip
gunzip gencode.v29.annotation.gtf.gz

# Sequin annotation
wget https://s3.amazonaws.com/sequins/annotations/RNAsequins.v2.2.gtf
# Checksum
md5sum RNAsequins.v2.2.gtf | grep f10f1da52825dccdf0cc23159ccf3735

# Combine annotations
cat gencode.v29.annotation.gtf RNAsequins.v2.2.gtf > GRCh38_spiked.gtf

# Remove individual files
rm RNAsequins.v2.2.gtf gencode.v29.annotation.gtf 

echo "Finished downloading human genome and sequin annotations"
