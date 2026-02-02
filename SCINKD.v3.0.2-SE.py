import os

#To run this pipeline on any machine running linux, run
#git clone https://github.com/DrPintoThe2nd/SCINKD.git
#mamba create -n scinkd meryl=1.4.1 snakemake=7.32.4 pigz r r-dplyr r-ggplot2 samtools minimap2 bbmap mosdepth --yes
#mamba activate scinkd
#mamba env export > SCINKD.v3.0.1a_environment.yml
#snakemake --use-conda --rerun-incomplete --nolock --cores 2 -j 1 -s SCINKD.v3.0.1a.py -np

configfile: "config.json"
R1 = config["R1_suffix"]
males = config["males"]
females = config["females"]
samples = males + females
genome = config["genome"]

rule all:
    input:
##setup links rule(s)
        expand("fastqs/{sample}.R1.fq.gz", sample=samples),
##meryl_count rule(s)
        expand("meryl_db/{sample}.meryl", sample=samples),
##meryl_intersect rule(s)
        "female_intersect.meryl",
        "male_intersect.meryl",
##meryl_diff rule(s)
        "female-specific.meryl/",
        "male-specific.meryl/",
##remove_noise
        "female-specific.1.meryl/",
        "male-specific.1.meryl/",
##meryl-lookup rule(s)
        expand("filtered_reads/{sample}.Fspec.R1.fastq.gz", sample=females),
        expand("filtered_reads/{sample}.Mspec.R1.fastq.gz", sample=males),
##calc_scinkd rule(s)
        expand("mapped_female/{sample}.bam", sample=females),
        expand("mapped_male/{sample}.bam", sample=males),
        expand("mapped_female/{sample}.bam.bai", sample=females),
        expand("mapped_male/{sample}.bam.bai", sample=males),
##coverage rule(s)
        expand("cov_F/{sample}.mosdepth.summary.txt", sample=females),
        expand("cov_M/{sample}.mosdepth.summary.txt", sample=males),
##results rule(s)

rule link_fastq_R1:
    input:
        "{sample}" + R1
    output:
        "fastqs/{sample}.R1.fq.gz"
    shell:
        """
        mkdir -p fastqs
        ln -sf $(realpath {input}) {output}
	    sleep 3
        """

rule meryl_count:
    input:
        r1 = "fastqs/{sample}.R1.fq.gz",
    output:
        out = directory("meryl_db/{sample}.meryl")
    params:
        threads = 24,
    shell:
        """
        mkdir -p meryl_db
        meryl count threads={params.threads} k=28 memory=24 {input.r1} output {output.out}
        sleep 3
        """

rule meryl_F_intersect:
    input:
        females = lambda wc: expand("meryl_db/{sample}.meryl", sample=females)
    output:
        directory("female_intersect.meryl")
    params:
        threads = 24,
    shell:
        """
        meryl intersect-min threads={params.threads} {input.females} output {output}
        """
rule meryl_M_intersect:
    input:
        males = lambda wc: expand("meryl_db/{sample}.meryl", sample=males)
    output:
        directory("male_intersect.meryl")
    params:
        threads = 24,
    shell:
        """
        meryl intersect-min threads={params.threads} {input.males} output {output}
        """

rule meryl_F_diff:
    input:
        F = "female_intersect.meryl",
        M = lambda wc: expand("meryl_db/{sample}.meryl", sample=males)
    output:
        directory("female-specific.meryl/"),
    shell:
        """
        meryl difference {input.F} {input.M} output {output}
        """

rule meryl_M_diff:
    input:
        M = "male_intersect.meryl",
        F = lambda wc: expand("meryl_db/{sample}.meryl", sample=females)
    output:
        directory("male-specific.meryl/"),
    shell:
        """
        meryl difference {input.M} {input.F} output {output}
        """

rule meryl_noise_F:
    input:
        "female-specific.meryl/",
    output:
        directory("female-specific.1.meryl/"),
    shell:
        """
        meryl greater-than 1 {input} output {output}
        """

rule meryl_noise_M:
    input:
        "male-specific.meryl/",
    output:
        directory("male-specific.1.meryl/"),
    shell:
        """
        meryl greater-than 1 {input} output {output}
        """

rule meryl_lookup_F1:
    input:
        r1 = "fastqs/{sample}.R1.fq.gz",
        mers = "female-specific.1.meryl/",
    output:
        "filtered_reads/{sample}.Fspec.R1.fastq.gz",
    shell:
        """
        mkdir -p filtered_reads/
        meryl-lookup -sequence {input.r1} -mers {input.mers} -include -output {output}
        """

rule meryl_lookup_F2:
    input:
        r2 = "fastqs/{sample}.R2.fq.gz",
        mers = "female-specific.1.meryl/",
    output:
        "filtered_reads/{sample}.Fspec.R2.fastq.gz",
    shell:
        """
        mkdir -p filtered_reads/
        meryl-lookup -sequence {input.r2} -mers {input.mers} -include -output {output}
        """

rule meryl_lookup_M1:
    input:
        r1 = "fastqs/{sample}.R1.fq.gz",
        mers = "male-specific.1.meryl/",
    output:
        "filtered_reads/{sample}.Mspec.R1.fastq.gz",
    shell:
        """
        mkdir -p filtered_reads/
        meryl-lookup -sequence {input.r1} -mers {input.mers} -include -output {output}
        """

##calculate number of kmers occuring in each haplotype
rule calc_scinkd_F:
    input:
        r1 = "filtered_reads/{sample}.Fspec.R1.fastq.gz",
    output:
        bam = "mapped_female/{sample}.bam",
        bai = "mapped_female/{sample}.bam.bai",
    params:
        threads = 12,
        genome = genome,
    shell:
        """
        mkdir -p mapped_female/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

rule calc_scinkd_M:
    input:
        r1 = "filtered_reads/{sample}.Mspec.R1.fastq.gz",
    output:
        bam = "mapped_male/{sample}.bam",
        bai = "mapped_male/{sample}.bam.bai",
    params:
        threads = 12,
        genome = genome,
    shell:
        """
        mkdir -p mapped_male/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

rule calc_cov_F:
    input:
        bam = "mapped_female/{sample}.bam",
    output:
        cov = "cov_F/{sample}.regions.bed.gz",
        sum = "cov_F/{sample}.mosdepth.summary.txt",
    params:
        threads = 4,
        prefix = lambda wc: "cov_F/" + wc.sample,
        window = 500000,
    shell:
        """
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """

rule calc_cov_M:
    input:
        bam = "mapped_male/{sample}.bam",
    output:
        cov = "cov_M/{sample}.regions.bed.gz",
        sum = "cov_M/{sample}.mosdepth.summary.txt",
    params:
        threads = 4,
        prefix = lambda wc: "cov_M/" + wc.sample,
        window = 500000,
    shell:
        """
        mkdir -p cov_M
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """
