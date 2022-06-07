# Assemble transcriptome and calculate gene/transcript abundance

## First follow the instructions here:
[Step by step guide on how to use my pipelines](https://carolinapb.github.io/2021-06-23-how-to-run-my-pipelines/)  
Click [here](https://github.com/CarolinaPB/snakemake-template/blob/master/Short%20introduction%20to%20Snakemake.pdf) for an introduction to Snakemake

## ABOUT

This is a pipeline that aligns raw RNA-seq reads (downloaded from SRA) to a genome and quantifies gene abundance.
Based on: _Pertea, M., Kim, D., Pertea, G. M., Leek, J. T., & Salzberg, S. L. (2016). Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie and Ballgown. Nature Protocols, 11(9), 1650–1667. https://doi.org/10.1038/nprot.2016.095_  
The final product is a table with gene abundances for every sample.

#### Tools used:
- [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) - download records from SRA
- [BBMAP (reformat.sh)](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/) - combine paired reads into one fastq
- [Hisat2](http://daehwankimlab.github.io/hisat2/manual/) - align RNA-seq to genome
- [Samtools](http://www.htslib.org/) - sort and index mapped reads
- [stringtie](https://ccb.jhu.edu/software/stringtie/#:~:text=StringTie%20is%20a%20fast%20and,variants%20for%20each%20gene%20locus.) - assemble alignments into transcripts

| ![DAG](https://github.com/CarolinaPB/get-expression-values/blob/master/workflow.png) |
|:--:|
|*Pipeline workflow* |

### Edit config.yaml with the paths to your files
```yaml
METALIST: /path/to/metalist.txt
ASSEMBLY:  /path/to/assembly.fa
ANNOTATION: /path/to/annotation.gff3
OUTDIR: /path/to/outdir
```

- METALIST - text file with sample details. Each row is a different sample/record. Mandatory columns:
    - Run - the most important - SRA code
    - AGE - age of individual. If not applicable, leave blank or add `NA`
    - tissue - tissue where the sample comes from. If not applicable, leave blank or add `NA`
- ASSEMBLY - decompressed fasta file 
- ANNOTATION - annotation gff3 file
- OUTDIR - directory where snakemake will run and where the results will be written to  
  If you want the results to be written to this directory (not to a new directory), open config.yaml and comment out `OUTDIR: /path/to/outdir`

## RESULTS
The most important results are:
- **<run_date>_files.txt** dated file with an overview of the files used to run the pipeline (for documentation purposes)
- **reads** directory that contains alignments
- **results** directory with final results
    - merged_transcriptome.gtf - non-redundant set of transcripts
    - results/{SRA}/{SRA}.gtf and results/{SRA}/gene_abundance.tab - re-estimated abundances for each sample
    - {SRA}/gene_abundance_{SRA}.txt - re-estimated abundances for each sample with extra informationa added - AGE and tissue
    - results/all_abundances.txt - final file with gene abundances for all samples

