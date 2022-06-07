# Assemble transcriptome and calculate gene/transcript abundance

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT

This is a pipeline that aligns raw RNA-seq reads (downloaded from SRA) to a genome and quantifies gene abundance.

#### Tools used:
- [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) - download records from SRA
- [BBMAP (reformat.sh)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - combine paired reads into one fastq
- [Hisat2](http://daehwankimlab.github.io/hisat2/manual/) - align RNA-seq to genome
- [Samtools](http://www.htslib.org/) - sort and index mapped reads
- [stringtie](https://ccb.jhu.edu/software/stringtie/#:~:text=StringTie%20is%20a%20fast%20and,variants%20for%20each%20gene%20locus.) - assemble alignments into transcripts

