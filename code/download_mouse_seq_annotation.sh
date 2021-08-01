#!/bin/bash

# Genome annotation
mkdir $data_dir/annotation_mouse
cd $data_dir/annotation_mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz
# Genome sequence
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz
# Checksum
md5sum gencode.vM21.annotation.gtf.gz | grep d08f66b2746d0ae66594fda6ea0c9939
md5sum GRCm38.primary_assembly.genome.fa.gz | grep e9b871a0de8039459595c2ed19545c2b

echo "Finished downloading mouse genome sequence and annotation"
