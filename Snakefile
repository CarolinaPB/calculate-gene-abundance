import pandas as pd
configfile: "config.yaml"

pipeline = "get-expression-values" # replace with your pipeline's name


include: "rules/create_file_log.smk"

METADATA = config["METALIST"]
ASSEMBLY = config["ASSEMBLY"]
ANNOTATION = config["ANNOTATION"]

RUNS = set()
for key, value in METADATA.items():
    key_val = pd.read_table(value, sep=",")
    run = set(key_val.Run)
    RUNS = RUNS | run

# print(RUNS)

def get_assembly_prefix(assembly):
    pref = assembly.rsplit('/', 1)[-1]
    prefix = pref.rsplit('.', 1)[0]
    return(prefix)

PREF = get_assembly_prefix(ASSEMBLY)


rule all:
    input:
        files_log,
        # expand("fastq/{SRA}.fastq", SRA=RUNS),
        # expand(PREF+"{ext}.ht2", ext= range(1,8))
        expand("{SRA}.gtf", SRA=RUNS)

# rule fastq_dump:
#     output:
#         "fastq/{SRA}.fastq"
#     message:
#         'Rule {rule} processing'
#     group:
#         'group'
#     params:
#         "{SRA}"
#     shell:
#         'fasterq-dump.2.11.0 -e 16 --outdir fastq {params}'

rule sam_dump:
    output:
        "reads/{SRA}.bam"
    message:
        'Rule {rule} processing'
    params:
        "{SRA}"
    group:
        'group'
    shell:
        'sam-dump {params} | samtools sort -@ 8 -O bam - > {output}'

rule hisat_build:
    input:
        ASSEMBLY
    output:
        expand(PREF+"{ext}.ht2", ext= range(1,8))
    message:
        'Rule {rule} processing'
    params:
        PREF
    shell:
        'hisat2-build {input} {params}'

rule stringtie:
    input:
        bam = rules.sam_dump.output,
        annotation = ANNOTATION,
        idx = rules.hisat_build.output
    output:
        "{SRA}.gtf"
    message:
        'Rule {rule} processing'
    group:
        'group'
    params:
        "{SRA}"
    shell:
        """
        stringtie {input.bam} \
        -G {input.annotation} \
        -o {params}.gtf \
        -p 16 \
        -A gene_abundance_{params}.tab \
        -B
        """
