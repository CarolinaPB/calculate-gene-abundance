configfile: "config.yaml"

import pandas as pd

pipeline = "get-expression-values" # replace with your pipeline's name

include: "rules/create_file_log.smk"

METADATA = config["METALIST"]
ASSEMBLY = config["ASSEMBLY"]
ANNOTATION = config["ANNOTATION"]

RUNS = set()
# for key, value in METADATA.items():
#     key_val = pd.read_table(value, sep=",")
#     run = set(key_val.Run)
#     RUNS = RUNS | run

metadata_table = pd.read_table(METADATA, sep=",")
run = set(metadata_table.Run)
RUNS = RUNS | run

# print(RUNS)

def get_assembly_prefix(assembly):
    pref = assembly.rsplit('/', 1)[-1]
    prefix = pref.rsplit('.', 1)[0]
    return(prefix)

PREF = get_assembly_prefix(ASSEMBLY)

localrules: create_file_log,fix_annotation_names

rule all:
    input:
        files_log,
        # expand("fastq/{SRA}.fastq", SRA=RUNS),
        # expand(PREF+"{ext}.ht2", ext= range(1,8))
        # expand("results/{SRA}/{SRA}.gtf", SRA=RUNS)
        "results/all_abundances.txt"

rule fastq_dump:
    output:
        "fastq/{SRA}.fastq"
    message:
        'Rule {rule} processing'
    group:
        'group'
    params:
        "{SRA}"
    shell:
        'fasterq-dump.2.11.0 -e 16 --outdir fastq {params}'

# rule sam_dump:
#     output:
#         "reads/{SRA}.bam"
#     message:
#         'Rule {rule} processing'
#     params:
#         "{SRA}"
#     group:
#         'group'
#     shell:
#         'sam-dump.2.11.0  {params} |  samtools view -bS - >  {output}'

rule hisat_build:
    input:
        ASSEMBLY
    output:
        expand(PREF+".{ext}.ht2", ext= range(1,9))
    message:
        'Rule {rule} processing'
    params:
        PREF
    shell:
        'hisat2-build {input} {params}'

rule hisat:
    input:
        rules.fastq_dump.output
    output:
        "reads/{SRA}.sam"
    message:
        'Rule {rule} processing'
    params:
        PREF
    group:
        'group'
    shell:
        # 'module load samtools && hisat2 -x {params} {input} | samtools sort -@ 8 - > {output}'
        'module load samtools && hisat2 -x {params} {input} > {output}'

rule sort_hisat:
    input:
        rules.hisat.output
    output:
        bam ="reads/{SRA}.bam",
        bai = "reads/{SRA}.bam.bai"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        '''
        module load samtools && samtools sort -@ 8 {input} > {output.bam}
        samtools index -@ 16 {output.bam}
        '''

# rule fix_annotation_names:
#     input:
#         ANNOTATION
#     output:
#         out = "data/annotation_fixed_names.gff3",
#         tmp1 = temp("tmp1"),
#         tmp2 = temp("tmp2")
#     message:
#         'Rule {rule} processing'
#     shell:
#         '''
#         cut -f9 {input} | awk -F'[=]' '{{print $2}}' | awk -F'[-]' '{{print $1}}' | awk -F'[;]' '{{print $1}}' > {output.tmp1}
#         cut -f -8 {input} > {output.tmp2}
#         paste {output.tmp2} {output.tmp1} > {output.out}
#         '''


rule stringtie:
    input:
        bam = rules.sort_hisat.output.bam,
        annotation=ANNOTATION,
        # annotation = rules.fix_annotation_names.output.out,
        idx = rules.hisat_build.output
    output:
        gtf = "results/{SRA}/{SRA}.gtf",
        abundance = "results/{SRA}/gene_abundance_{SRA}.tab"
    message:
        'Rule {rule} processing'
    group:
        'group'
    # params:
    #     "{SRA}"
    shell:
        """
        stringtie {input.bam} \
        -G {input.annotation} \
        -o {output.gtf} \
        -p 16 \
        -A {output.abundance} \
        -B
        """

rule concat_res:
    input:
        expand("results/{SRA}/gene_abundance_{SRA}.tab", SRA=RUNS)
    output:
        "results/all_abundances.txt"
    message:
        'Rule {rule} processing'
    group:
        'group'
    shell:
        'cat {input} > {output}'

