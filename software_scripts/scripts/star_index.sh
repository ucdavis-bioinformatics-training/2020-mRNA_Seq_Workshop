#!/bin/bash

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz
FASTA="../GRCh38.primary_assembly.genome.fa"

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
gunzip gencode.v34.primary_assembly.annotation.gtf.gz
GTF="../gencode.v34.primary_assembly.annotation.gtf"

mkdir star.overlap100.gencode.v34
cd star.overlap100.gencode.v34

module load star

call="STAR
     --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir . \
     --sjdbOverhang 100 \
     --sjdbGTFfile ${GTF} \
     --genomeFastaFiles ${FASTA}"

echo $call
eval $call

end=`date +%s`
runtime=$((end-start))
echo $runtime
