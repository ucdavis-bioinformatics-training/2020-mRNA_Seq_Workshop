#!/bin/bash

start=`date +%s`
echo $HOSTNAME

outpath="References"
mkdir -p ${outpath}

cd ${outpath}
wget ftp://ftp.ensembl.org/pub/release-100/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz
gunzip Mus_musculus.GRCm38.dna.toplevel.fa.gz
FASTA="../Mus_musculus.GRCm38.dna.toplevel.fa"

wget ftp://ftp.ensembl.org/pub/release-100/gtf/mus_musculus/Mus_musculus.GRCm38.100.gtf.gz
gunzip Mus_musculus.GRCm38.100.gtf.gz
GTF="../Mus_musculus.GRCm38.100.gtf"

mkdir star_2.7.3a_index_GRCm38.p6
cd star_2.7.3a_index_GRCm38.p6


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
