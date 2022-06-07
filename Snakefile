#################################
# author: Carolina Pita Barros  #
# carolina.pitabarros@wur.nl    #
# date: July 2021            #
#################################

configfile: "config.yaml"

from snakemake.utils import makedirs
import pandas as pd
from pathlib import Path

pipeline = "get-expression-values"  


include: "rules/create_file_log.smk"


METADATA = config["METALIST"]
ASSEMBLY = config["ASSEMBLY"]
ANNOTATION = config["ANNOTATION"]

if "OUTDIR" in config:
    workdir: config["OUTDIR"]
makedirs("logs_slurm")

# Get SRA codes of data to download
metadata_table = pd.read_table(METADATA, sep=",")
RUNS = set(metadata_table.Run)

PREF = Path(ASSEMBLY).stem


localrules:
    create_file_log,
    add_info_to_stringtie,
    create_list,


rule all:
    input:
        files_log,
        expand("results/{SRA}/gene_abundance_{SRA}.txt", SRA=RUNS),
        "results/all_abundances.txt",

rule fastq_dump:    
    """
    Download data from SRA
    """
    output:
        temp("fastq/{SRA}.fastq"),
    message:
        "Rule {rule} processing"
    group:
        "group_{SRA}"
    params:
        "{SRA}",
    shell:
        "fasterq-dump.2.11.0 -e 16 --outdir fastq {params}"

rule combine_fastq:
    """
    If a sample has forward and reverse reads, merge both into one unique fastq
    """
    input:
        in1 = 'fastq/{SRA}_1.fastq',
        in2 = 'fastq/{SRA}_2.fastq'
    output:
        temp('fastq/{SRA}.fastq')
    message:
        'Rule {rule} processing'
    shell:
        'reformat.sh in1={input.in1} in2={input.in2} out={output}'

rule hisat_build:
    """
    Index genome
    """
    input:
        ASSEMBLY,
    output:
        expand(f"{PREF}.{{ext}}.ht2", ext=range(1, 9)),
    message:
        "Rule {rule} processing"
    params:
        PREF,
    shell:
        "hisat2-build {input} {params}"


rule hisat:
    """
    Align RNA-seq reads to genome
    """
    input:
        fq=rules.fastq_dump.output,
        idx=rules.hisat_build.output,
    output:
        temp("reads/{SRA}.sam"),
    message:
        "Rule {rule} processing"
    params:
        PREF,
    group:
        "group_{SRA}"
    shell:
        "hisat2 -x {params} {input.fq} > {output}"


rule sort_hisat:
    input:
        rules.hisat.output,
    output:
        bam="reads/{SRA}.bam",
        bai="reads/{SRA}.bam.bai",
    message:
        "Rule {rule} processing"
    group:
        "group_{SRA}"
    shell:
        """
        module load samtools
        samtools sort -@ 16 {input} > {output.bam}
        samtools index -@ 16 {output.bam}
        """


rule stringtie:
    """
    Assemble RNA-seq alignments into potential transcripts
    """
    input:
        bam=rules.sort_hisat.output.bam,
        annotation=ANNOTATION,
    output:
        gtf=temp("before_merge/{SRA}/{SRA}.gtf"),
        abundance=temp("before_merge/{SRA}/gene_abundance.tab"),
    message:
        "Rule {rule} processing"
    group:
        "group_{SRA}"
    params:
        "{SRA}",
    shell:
        """
        stringtie {input.bam} \
        -G {input.annotation} \
        -l {params} \
        -o {output.gtf} \
        -p 16 \
        -A {output.abundance} \
        -B \
        -e 
        """

rule create_list:
    """
    Create list with gtf files to merge
    """
    input:
        expand("before_merge/{SRA}/{SRA}.gtf", SRA=RUNS),
    output:
        temp("list_to_merge.txt"),
    message:
        "Rule {rule} processing"
    run:
        with open(output[0], "w") as out:
            for el in input:
                out.write(el)
                out.write("\n")


rule stringtie_merge:
    """
    Create non-redundant set of trasncripts from list of transcripts. 
    Creates one set of transcripts for all RNA-seq samples.
    """
    input:
        list_to_merge=rules.create_list.output,
        annotation=ANNOTATION,
    output:
        "results/merged_transcriptome.gtf",
    message:
        "Rule {rule} processing"
    shell:
        """
        stringtie --merge \
        -G {input.annotation} \
        -o {output} {input.list_to_merge}
        """


rule stringtie_after_merge:
    """
    Re-estimate transcript abundances for each sample
    """
    input:
        bam=rules.sort_hisat.output.bam,
        annotation=rules.stringtie_merge.output,
    output:
        gtf="results/{SRA}/{SRA}.gtf",
        abundance="results/{SRA}/gene_abundance.tab",
    message:
        "Rule {rule} processing"
    group:
        "group_{SRA}"
    params:
        "{SRA}",
    shell:
        """
        stringtie {input.bam} \
        -G {input.annotation} \
        -l {params} \
        -o {output.gtf} \
        -p 16 \
        -A {output.abundance} \
        -B \
        -e 
        """


rule add_info_to_stringtie:
    """
    Add run, age and tissue data to stringtie results
    """
    input:
        rules.stringtie_after_merge.output.abundance,
        METADATA
    output:
        "results/{SRA}/gene_abundance_{SRA}.txt",
    message:
        "Rule {rule} processing"
    group:
        "group_{SRA}"
    params:
        "{SRA}",
    run:
        abundance = pd.read_table(input[0])
        metadata = pd.read_table(input[1], sep=",")
        details = metadata[metadata["Run"] == params[0]]
        details_subset = details[["Run", "AGE", "tissue"]]
        details_subset.reset_index(inplace=True)
        combined_tables = pd.concat([abundance, details_subset], axis=1)
        combined_tables["Run"] = combined_tables.iloc[0]["Run"]
        combined_tables["AGE"] = combined_tables.iloc[0]["AGE"]
        combined_tables["tissue"] = combined_tables.iloc[0]["tissue"]
        combined_tables.to_csv(output[0], sep="\t", index=False)


rule concat_res:
    """
    Combine abundance results from each sample into one final file
    """
    input:
        expand("results/{SRA}/gene_abundance_{SRA}.txt", SRA=RUNS),
    output:
        "results/all_abundances.txt",
    message:
        "Rule {rule} processing"
    shell:
        "cat {input} > {output}"
