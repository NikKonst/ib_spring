import pandas as pd
import os

configfile: 'workflow/config.yaml'

SAMPLES_INFO = pd.read_csv(config['data_info'], sep=';')
SAMPLES_INFO.Target = SAMPLES_INFO.Target.fillna('')
SAMPLES_INFO['Samples'] = SAMPLES_INFO.GSM + '_' + SAMPLES_INFO.Cell + '_' + SAMPLES_INFO.Method
SAMPLES_INFO = SAMPLES_INFO.set_index(['Samples'], drop=False)

SAMPLES_INFO_INDEX = SAMPLES_INFO.index.tolist()
ind_n = 1
for i, b in enumerate(SAMPLES_INFO.index.duplicated()):
    if b:
        SAMPLES_INFO_INDEX[i] += '_' + str(ind_n)
        ind_n += 1
    else:
        ind_n = 1

SAMPLES_INFO.index = SAMPLES_INFO_INDEX

SAMPLES_INFO['File'] = SAMPLES_INFO.SRX + '/' + SAMPLES_INFO.SRR + '.fastq'

PEAK_CALLING_RES = ["_peaks.xls", "_peaks.narrowPeak", "_summits.bed"]

rule all:
    input:
        # expand(f"data/reads/{{file}}", file=SAMPLES_INFO.File)
        'qc/multiqc/reads.html',
        'qc/multiqc/bams.html',
        expand(f"bigwigs/{{sample}}_{config['genome']}.bw", sample=SAMPLES_INFO.Samples),
        # expand(f"callpeak/{{sample}}_{config['genome']}{{ext}}", sample=SAMPLES_INFO.Samples, ext=PEAK_CALLING_RES)

rule get_fastq:
    output:
        "data/reads/{srx}/{accession}.fastq"
    params:
        tmp_dir = config['tmp_dir']
    log:
        'logs/fastq_dump/{srx}/{accession}.log'
    threads: 8
    run:
        out_dir = os.path.dirname(output[0])
        shell("parallel-fastq-dump --threads {threads} "
              f"--outdir {out_dir} --sra-id {{wildcards.accession}}")

rule fastqc:
    input:
        lambda wildcards: f"data/reads/{SAMPLES_INFO.loc[wildcards, 'File'][0]}"
    output:
        html='qc/fastqc/{sample}.html',
        zip='qc/fastqc/{sample}_fastqc.zip'
    params: '--quiet'
    log:
        'logs/fastqc/{sample}.log'
    threads: 8
    wrapper:
        '0.73.0/bio/fastqc'

rule multiqc:
    input:
        expand('qc/fastqc/{sample}_fastqc.zip', sample=SAMPLES_INFO.Samples)
    output:
        'qc/multiqc/reads.html'
    log:
        'logs/multiqc/multiqc.log'
    wrapper:
        '0.73.0/bio/multiqc'

rule download_reference_genome:
    output:
        f"indexes/{config['genome']}/{config['genome']}.fa.gz"
    shell:
        f"wget http://hgdownload.soe.ucsc.edu/goldenPath/{config['genome']}/bigZips/{config['genome']}.fa.gz -O {{output}}"

rule unzip_reference_genome:
    input:
        f"indexes/{config['genome']}/{config['genome']}.fa.gz"
    output:
        f"indexes/{config['genome']}/{config['genome']}.fa"
    shell:
        "gzip -d -c {input} > {output}"

rule bowtie2_build:
    input:
        reference=f"indexes/{config['genome']}/{config['genome']}.fa"
    output:
        multiext(
            f"indexes/{config['genome']}/{config['genome']}",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    conda:
        'envs/bowtie2.yaml'
    log:
        "logs/bowtie2_build/build.log"
    threads: 8
    script:
        'scripts/bowtie2-build.py'

rule bowtie2_align:
    input:
        sample=lambda wildcards: [f"data/reads/{SAMPLES_INFO.loc[wildcards, 'File'][0]}"],
        index_files=multiext(
            f"indexes/{config['genome']}/{config['genome']}",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    output:
        f"bams/{{sample}}_{config['genome']}.bam"
    conda:
        'envs/bowtie2.yaml'
    log:
        "logs/bowtie2_align/{sample}.log"
    params:
        index=f"indexes/{config['genome']}/{config['genome']}"
    threads: 8
    script:
        'scripts/bowtie2-align.py'

rule multiqc_bowtie2:
    input:
        expand(f"logs/bowtie2_align/{{sample}}.log", sample=SAMPLES_INFO.Samples)
    output:
        'qc/multiqc/bams.html'
    log:
        'logs/multiqc/multiqc_bowtie2.log'
    wrapper:
        '0.73.0/bio/multiqc'

rule samtools_sort_bam:
    input:
        f"bams/{{sample}}_{config['genome']}.bam"
    output:
        f"bams/{{sample}}_{config['genome']}.sorted.bam"
    params:
        extra = "-m 4G",
        tmp_dir = config['tmp_dir']
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.73.0/bio/samtools/sort"

rule samtools_index_bam:
    input:
        f"bams/{{sample}}_{config['genome']}.sorted.bam"
    output:
        f"bams/{{sample}}_{config['genome']}.sorted.bam.bai"
    wrapper:
        "0.73.0/bio/samtools/index"

rule bam_coverage:
    input:
        bams=f"bams/{{sample}}_{config['genome']}.sorted.bam",
        bam_indexes=f"bams/{{sample}}_{config['genome']}.sorted.bam.bai"
    output:
        f"bigwigs/{{sample}}_{config['genome']}.bw"
    conda:
        config['conda_env']
    shell:
        f"bamCoverage --bam {{input.bams}} -o {{output}}"

# rule callpeak:
#     input:
#         treatment=f"bams/{{sample}}_{config['genome']}.sorted.bam",   # required: treatment sample(s)
#     output:
#         multiext(f"callpeak/{{sample}}_{config['genome']}",
#                  "_peaks.xls",   ### required
#                  ### optional output files
#                  "_peaks.narrowPeak",
#                  "_summits.bed"
#                  )
#     log:
#         "logs/macs2/callpeak.{sample}.log"
#     params:
#         "-f BAM -g hs --nomodel"
#     wrapper:
#         "0.73.0/bio/macs2/callpeak"
