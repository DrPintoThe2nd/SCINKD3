import os

#To run this pipeline on any machine running linux, run
#git clone https://github.com/DrPintoThe2nd/SCINKD.git
#mamba create -n scinkd meryl=1.4.1 snakemake=7.32.4 pigz r r-dplyr r-ggplot2 samtools minimap2 bbmap mosdepth --yes
#mamba activate scinkd
#mamba env export > SCINKD.v3.0.1a_environment.yml
#snakemake --use-conda --rerun-incomplete --nolock --cores 2 -j 1 -s SCINKD.v3.0.2.py -np

configfile: "SCINKD3/SCINKD.3.1.0_config.json"
R1 = config["R1_suffix"]
R2 = config["R2_suffix"]
males = config["males"]
females = config["females"]
samples = males + females
genome = config["genome"]

threads = config["per_job_threads"]
meryl_mem = config["meryl_mem"]
repair_mem = config["repair_mem"]
window = config["depth_window"]
ChrNum = config["ChrNum"]

rule all:
    input:
##setup links rule(s)
        expand("fastqs/{sample}.R1.fq.gz", sample=samples),
        expand("fastqs/{sample}.R2.fq.gz", sample=samples),
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
        expand("filtered_reads/{sample}.Fspec.R2.fastq.gz", sample=females),
        expand("filtered_reads/{sample}.Mspec.R1.fastq.gz", sample=males),
        expand("filtered_reads/{sample}.Mspec.R2.fastq.gz", sample=males),
##repair rule(s)
        expand("repaired_reads/{sample}.Fspec.R1.fq.gz", sample=females),
        expand("repaired_reads/{sample}.Fspec.R2.fq.gz", sample=females),
        expand("repaired_reads/{sample}.Mspec.R1.fq.gz", sample=males),
        expand("repaired_reads/{sample}.Mspec.R2.fq.gz", sample=males),
##cat rule(s)
        "repaired_reads/Fspec.R1.fastq.gz",
        "repaired_reads/Fspec.R2.fastq.gz",
        "repaired_reads/Mspec.R1.fastq.gz",
        "repaired_reads/Mspec.R2.fastq.gz",
##mapping_scinkd rule(s)
        expand("mapped_female/{sample}.bam", sample=females),
        expand("mapped_male/{sample}.bam", sample=males),
        expand("mapped_female/{sample}.bam.bai", sample=females),
        expand("mapped_male/{sample}.bam.bai", sample=males),
        "mapped_female/Fspec.bam.bai",
        "mapped_male/Mspec.bam.bai",
##coverage rule(s)
        expand("cov_F/F_{sample}.mosdepth.summary.txt", sample=females),
        expand("cov_M/M_{sample}.mosdepth.summary.txt", sample=males),
        "Fspec_total.regions.bed.gz",
        "Fspec_total.mosdepth.summary.txt",
        "Mspec_total.regions.bed.gz",
        "Mspec_total.mosdepth.summary.txt",
##individual results rule(s)
        expand("cov_F/F_{sample}.regions.manhattan.png", sample=females),
        expand("cov_F/F_{sample}.regions.manhattan.pdf", sample=females),
        expand("cov_M/M_{sample}.regions.manhattan.png", sample=males),
        expand("cov_M/M_{sample}.regions.manhattan.pdf", sample=males),
##combined result(s)
        "Fspec_total.regions.manhattan.png",
        "Mspec_total.regions.manhattan.pdf",
#        expand("{genome}.dotplot.png", genome=genome),
#        expand("{genome}.dotplot.pdf", genome=genome),

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

rule link_fastq_R2:
    input:
        "{sample}" + R2
    output:
        "fastqs/{sample}.R2.fq.gz"
    shell:
        """
        mkdir -p fastqs
        ln -sf $(realpath {input}) {output}
	    sleep 2
        """

rule meryl_count:
    input:
        r1 = "fastqs/{sample}.R1.fq.gz",
        r2 = "fastqs/{sample}.R2.fq.gz"
    output:
        out = directory("meryl_db/{sample}.meryl")
    params:
        threads = threads,
        memory = meryl_mem,
    shell:
        """
        mkdir -p meryl_db
        meryl count threads={params.threads} k=28 memory={params.memory} {input.r1} {input.r2} output {output.out}
        sleep 2
        """

rule meryl_F_intersect:
    input:
        females = lambda wc: expand("meryl_db/{sample}.meryl", sample=females)
    output:
        directory("female_intersect.meryl")
    params:
        threads = threads,
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
        threads = threads,
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

rule meryl_lookup_M2:
    input:
        r2 = "fastqs/{sample}.R2.fq.gz",
        mers = "male-specific.1.meryl/",
    output:
        "filtered_reads/{sample}.Mspec.R2.fastq.gz",
    shell:
        """
        mkdir -p filtered_reads/
        meryl-lookup -sequence {input.r2} -mers {input.mers} -include -output {output}
        """

rule repair_F_reads:
    input:
        r1 = "filtered_reads/{sample}.Fspec.R1.fastq.gz",
        r2 = "filtered_reads/{sample}.Fspec.R2.fastq.gz",
    output:
        r1 = "repaired_reads/{sample}.Fspec.R1.fq.gz",
        r2 = "repaired_reads/{sample}.Fspec.R2.fq.gz",
        S  = "repaired_reads/{sample}.Fspec.S.fastq.gz",
    params:
        memory = repair_mem,
    shell:
        """
        mkdir -p repaired_reads/
        repair.sh -Xmx{params.memory}g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.S} -da showspeed=t overwrite=t 2>/dev/null
        """

rule repair_M_reads:
    input:
        r1 = "filtered_reads/{sample}.Mspec.R1.fastq.gz",
        r2 = "filtered_reads/{sample}.Mspec.R2.fastq.gz",
    output:
        r1 = "repaired_reads/{sample}.Mspec.R1.fq.gz",
        r2 = "repaired_reads/{sample}.Mspec.R2.fq.gz",
        S  = "repaired_reads/{sample}.Mspec.S.fastq.gz",
    params:
        memory = repair_mem,
    shell:
        """
        mkdir -p repaired_reads/
        repair.sh -Xmx{params.memory}g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.S} -da showspeed=t overwrite=t 2>/dev/null
        """

##map unique kmers from in each individual
rule calc_scinkd_F:
    input:
        r1 = "repaired_reads/{sample}.Fspec.R1.fq.gz",
        r2 = "repaired_reads/{sample}.Fspec.R2.fq.gz",
    output:
        bam = "mapped_female/{sample}.bam",
        bai = "mapped_female/{sample}.bam.bai",
    params:
        threads = threads,
        genome = genome,
    shell:
        """
        mkdir -p mapped_female/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

rule calc_scinkd_M:
    input:
        r1 = "repaired_reads/{sample}.Mspec.R1.fq.gz",
        r2 = "repaired_reads/{sample}.Mspec.R2.fq.gz",
    output:
        bam = "mapped_male/{sample}.bam",
        bai = "mapped_male/{sample}.bam.bai",
    params:
        threads = threads,
        genome = genome,
    shell:
        """
        mkdir -p mapped_male/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

#calculate coverage/depth for each individual
rule calc_cov_F:
    input:
        bam = "mapped_female/{sample}.bam",
    output:
        cov = "cov_F/F_{sample}.regions.bed.gz",
        sum = "cov_F/F_{sample}.mosdepth.summary.txt",
    params:
        threads = threads,
        prefix = lambda wc: "cov_F/F_" + wc.sample,
        window = window,
    shell:
        """
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """

rule calc_cov_M:
    input:
        bam = "mapped_male/{sample}.bam",
    output:
        cov = "cov_M/M_{sample}.regions.bed.gz",
        sum = "cov_M/M_{sample}.mosdepth.summary.txt",
    params:
        threads = threads,
        prefix = lambda wc: "cov_M/M_" + wc.sample,
        window = window,
    shell:
        """
        mkdir -p cov_M
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """

rule plot_scinkd_F:
    input:
        cov = "cov_F/F_{sample}.regions.bed.gz",
    output:
        png = "cov_F/F_{sample}.regions.manhattan.png",
        pdf = "cov_F/F_{sample}.regions.manhattan.pdf",
    params:
        ChrNum = ChrNum,
    shell:
        """
        Rscript manhattan_plot_indiv.R {input.cov} {params.ChrNum}
        """

rule plot_scinkd_M:
    input:
        cov = "cov_M/M_{sample}.regions.bed.gz",
    output:
        png = "cov_M/M_{sample}.regions.manhattan.png",
        pdf = "cov_M/M_{sample}.regions.manhattan.pdf",
    params:
        ChrNum = ChrNum,
    shell:
        """
        Rscript manhattan_plot_indiv.R {input.cov} {params.ChrNum}
        """

##combine data from each individual to summarize and amplify signal
##first combine reads
##calculate number of kmers occuring in each haplotype
rule cat_Fspec:
    input:
        r1 = expand("repaired_reads/{sample}.Fspec.R1.fq.gz", sample=females),
        r2 = expand("repaired_reads/{sample}.Fspec.R2.fq.gz", sample=females),
    output:
        cat1 = "repaired_reads/Fspec.R1.fastq.gz",
        cat2 = "repaired_reads/Fspec.R2.fastq.gz",
    shell:
        """
        cat {input.r1} > {output.cat1}
        cat {input.r2} > {output.cat2}        
        """

rule cat_Mspec:
    input:
        r1 = expand("repaired_reads/{sample}.Mspec.R1.fq.gz", sample=males),
        r2 = expand("repaired_reads/{sample}.Mspec.R2.fq.gz", sample=males),
    output:
        cat1 = "repaired_reads/Mspec.R1.fastq.gz",
        cat2 = "repaired_reads/Mspec.R2.fastq.gz",
    shell:
        """
        cat {input.r1} > {output.cat1}
        cat {input.r2} > {output.cat2}        
        """

rule mapcat_scinkd_F:
    input:
        r1 = "repaired_reads/Fspec.R1.fastq.gz",
        r2 = "repaired_reads/Fspec.R2.fastq.gz",
    output:
        bam = "mapped_female/Fspec.bam",
        bai = "mapped_female/Fspec.bam.bai",
    params:
        threads = threads,
        genome = genome,
    shell:
        """
        mkdir -p mapped_female/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

rule mapcat_scinkd_M:
    input:
        r1 = "repaired_reads/Mspec.R1.fastq.gz",
        r2 = "repaired_reads/Mspec.R2.fastq.gz",
    output:
        bam = "mapped_male/Mspec.bam",
        bai = "mapped_male/Mspec.bam.bai",
    params:
        threads = threads,
        genome = genome,
    shell:
        """
        mkdir -p mapped_male/
        minimap2 -ax sr -t{params.threads} {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{params.threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
        """

#calculate coverage/depth for combined throughput
rule cov_Fspec:
    input:
        bam = "mapped_female/Fspec.bam",
    output:
        cov = "Fspec_total.regions.bed.gz",
        sum = "Fspec_total.mosdepth.summary.txt",
    params:
        threads = threads,
        prefix = "Fspec_total",
        window = window,
    shell:
        """
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """

rule cov_Mspec:
    input:
        bam = "mapped_male/Mspec.bam",
    output:
        cov = "Mspec_total.regions.bed.gz",
        sum = "Mspec_total.mosdepth.summary.txt",
    params:
        threads = threads,
        prefix = "Mspec_total",
        window = window,
    shell:
        """
        MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{params.threads} {params.prefix} {input.bam}
        """

rule plot_Fspec:
    input:
        cov = "Fspec_total.regions.bed.gz",
    output:
        png = "Fspec_total.regions.manhattan.png",
        pdf = "Fspec_total.regions.manhattan.pdf",
    params:
        ChrNum = ChrNum,
    shell:
        """
        Rscript manhattan_plot_indiv.R {input.cov} {params.ChrNum}
        """

rule plot_Mspec:
    input:
        cov = "Mspec_total.regions.bed.gz",
    output:
        png = "Mspec_total.regions.manhattan.png",
        pdf = "Mspec_total.regions.manhattan.pdf",
    params:
        ChrNum = ChrNum,
    shell:
        """
        Rscript manhattan_plot_indiv.R {input.cov} {params.ChrNum}
        """

#rule plot_scinkd:
#    input:
#        F = "Fspec_total.regions.bed.gz",
#        M = "Mspec_total.regions.bed.gz",
#    output:
#        png = "{genome}.dotplot.png",
#        pdf = "{genome}.dotplot.pdf",
#    params:
#        ChrNum = ChrNum,
#        genome = genome,
#    shell:
#        """
#        Rscript chrom_dotplot.R {params.genome} {input.F} {input.M} {params.ChrNum}
#        """
#
