from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)
indexbase = snakemake.output[0].replace(".1.bt2", "")
shell(
    "bowtie2-build --threads {snakemake.threads} "
    "{snakemake.input.reference} {indexbase}"
)