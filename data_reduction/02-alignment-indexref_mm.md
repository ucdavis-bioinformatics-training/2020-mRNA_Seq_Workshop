
# Indexing a Reference sequence and annotation


This document assumes [preproc htstream](./01-preproc_htstream_mm.md) has been completed.

**IF** for some reason it didn't finish, is corrupted or you missed the session, you can copy over a completed copy

```bash
cp -r /share/biocore/workshops/2020_mRNAseq_July/HTS_testing /share/workshop/mrnaseq_workshop/$USER/rnaseq_example/.
cp -r /share/biocore/workshops/2020_mRNAseq_July/01-HTS_Preproc /share/workshop/mrnaseq_workshop/$USER/rnaseq_example/.
```

# Indexing the reference.

Before we can align reads to a reference genome/transcriptome, its typical to first build an index. Indexing a genome/transcriptome is similar to indexing a book. If you want to know on which page a certain word appears or a chapter begins (chromosome section a kmer appears in our case), it is much more efficient/faster to look it up in a pre-built index than going through every page of the book until you find it. Same goes for alignments. Indices allow the aligner to narrow down the potential locations of a query sequence within the genome, saving both time and memory.

1. First lets make sure we are where we are supposed to be and that the References directory is available.

    ```bash
    cd /share/workshop/mrnaseq_workshop/$USER/rnaseq_example
    mkdir slurmout
    ```

1. To align our data we will need the genome (fasta) and annotation (gtf) for mouse. There are many places to find them, but we are going to get them from the [Ensembl](https://uswest.ensembl.org/Mus_musculus/Info/Index).

    We need to first get the url for the genome and annotation gtf. For RNAseq we want to use the primary genome chromosomes and basic gene annotation. At the time of this workshop the current version of the Ensembl mouse genome is GRCm38.p6. You will want to update the scripts to use the current version.

    <img src="alignment_figures/ensembl1.png" alt="index_figure1" width="80%" style="border:5px solid #ADD8E6;"/>

1. We are going to use an aligner called ['STAR'](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) to align the data. Lets take a look at the help docs for star:

    ```bash
    module load star
    STAR -h
    ```

1. First we need to index the genome for STAR. Lets pull down a slurm script to index the Ensembl version of the mouse genome.

    ```bash
    wget https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2020-mRNA_Seq_Workshop/master/software_scripts/scripts/star_index.slurm
    less star_index.slurm
    ```

    <pre class="prettyprint"><code class="language-py" style="background-color:333333">
        #!/bin/bash
        #SBATCH --job-name=star_index # Job name
        #SBATCH --nodes=1
        #SBATCH --ntasks=8
        #SBATCH --time=120
        #SBATCH --mem=40000 # Memory pool for all cores (see also --mem-per-cpu)
        #SBATCH --partition=production
        #SBATCH --reservation=mrnaseq_workshop
        #SBATCH --account=mrnaseq_workshop
        #SBATCH --output=slurmout/star-index_%A.out # File to which STDOUT will be written
        #SBATCH --error=slurmout/star-index_%A.err # File to which STDERR will be written
        #SBATCH --mail-type=ALL
        #SBATCH --mail-user=myemail@email.com

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
    </code></pre>

{:start="1"}
1. The script uses wget to download the fasta and GTF files from Ensembl using the links you found earlier.
1. Uncompresses them using gunzip.
1. Creates the star index directory [star_2.7.3a_index_GRCm38.p6].
1. Change directory into the new star index directory. We run the star indexing command from inside the directory, for some reason star fails if you try to run it outside this directory.
1. Run star in mode genomeGenerate.

When you are done, type "q" to exit.

{:start="5"}
1. Run star indexing when ready.

    ```bash
    sbatch star_index.slurm
    ```

    This step will take a couple hours. You can look at the [STAR documentation](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) while you wait. All of the output files will be written to the star index directory **star_2.7.3a_index_GRCm38.p6**.

    **IF** For the sake of time, or for some reason it didn't finish, is corrupted, or you missed the session, you can **link** over a completed copy.

    ```bash
    ln -s /share/biocore/workshops/2020_mRNAseq_July/References/star_2.7.3a_index_GRCm38.p6 /share/workshop/mrnaseq_workshop/$USER/rnaseq_example/References/.
    ```
