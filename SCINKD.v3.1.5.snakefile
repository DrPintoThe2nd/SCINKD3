import os

#To run this pipeline on any machine running linux, run
#git clone https://github.com/DrPintoThe2nd/SCINKD3.git
#mamba create -n scinkd meryl=1.4.1 snakemake=7.32.4 pigz r r-dplyr r-ggplot2 samtools minimap2 bbmap mosdepth --yes
#mamba activate scinkd
#mamba env export > SCINKD.v3.0.1a_environment.yml
#snakemake --use-conda --rerun-incomplete --nolock --cores 2 -j 1 -s SCINKD.v3.0.2.py -np

configfile: "SCINKD3/config_SCINKD.v3.1.5.json"
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
		expand("filtered_reads/females/{sample}.Fspec.R1.fastq.gz", sample=females),
		expand("filtered_reads/females/{sample}.Fspec.R2.fastq.gz", sample=females),
		expand("filtered_reads/males/{sample}.Mspec.R1.fastq.gz", sample=males),
		expand("filtered_reads/males/{sample}.Mspec.R2.fastq.gz", sample=males),
##repair rule(s)
		expand("repaired_reads/females/{sample}.Fspec.R1.fq.gz", sample=females),
		expand("repaired_reads/females/{sample}.Fspec.R2.fq.gz", sample=females),
		expand("repaired_reads/males/{sample}.Mspec.R1.fq.gz", sample=males),
		expand("repaired_reads/males/{sample}.Mspec.R2.fq.gz", sample=males),
##cat rule(s)
		"repaired_reads/Fspec.R1.fastq.gz",
		"repaired_reads/Fspec.R2.fastq.gz",
		"repaired_reads/Mspec.R1.fastq.gz",
		"repaired_reads/Mspec.R2.fastq.gz",
##mapping_scinkd rule(s)
		expand("mapped/females/{sample}.bam", sample=females),
		expand("mapped/females/{sample}.bam.bai", sample=females),
		expand("mapped/males/{sample}.bam", sample=males),
		expand("mapped/males/{sample}.bam.bai", sample=males),
		"mapped/Fspec.bam",
		"mapped/Fspec.bam.bai",
		"mapped/Mspec.bam",
		"mapped/Mspec.bam.bai",
##coverage rule(s)
		expand("coverage/females/{sample}.mosdepth.summary.txt", sample=females),
		expand("coverage/males/{sample}.mosdepth.summary.txt", sample=males),
		"coverage/Fspec_total.regions.bed.gz",
		"coverage/Fspec_total.mosdepth.summary.txt",
		"coverage/Mspec_total.regions.bed.gz",
		"coverage/Mspec_total.mosdepth.summary.txt",
##individual results rule(s)
		expand("coverage/females/{sample}.regions.manhattan.png", sample=females),
		expand("coverage/females/{sample}.regions.manhattan.pdf", sample=females),
		expand("coverage/males/{sample}.regions.manhattan.png", sample=males),
		expand("coverage/males/{sample}.regions.manhattan.pdf", sample=males),
##combined result(s)
		"coverage/Fspec_total.regions.manhattan.png",
		"coverage/Mspec_total.regions.manhattan.pdf",
		"SCINKD3.dotplot.png",
		"SCINKD3.dotplot.pdf",
		"SCINKD3.females.png",
		"SCINKD3.females.pdf",
		"SCINKD3.males.png",
		"SCINKD3.males.pdf",

#homogenize read names for analysis
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

#count kmers in each sample
rule meryl_count:
	input:
		r1 = "fastqs/{sample}.R1.fq.gz",
		r2 = "fastqs/{sample}.R2.fq.gz"
	output:
		out = directory("meryl_db/{sample}.meryl")
	resources:
		mem_gb = meryl_mem,
	threads: threads,
	shell:
		"""
		mkdir -p meryl_db
		meryl count threads={threads} k=28 memory={resources.mem_gb} {input.r1} {input.r2} output {output.out}
		sleep 2
		"""

#intersect sex-specific kmers
rule meryl_intersect_F:
	input:
		females = lambda wc: expand("meryl_db/{sample}.meryl", sample=females)
	output:
		directory("female_intersect.meryl")
	threads: threads,
	shell:
		"""
		meryl intersect-min threads={threads} {input.females} output {output}
		"""
rule meryl_intersect_M:
	input:
		males = lambda wc: expand("meryl_db/{sample}.meryl", sample=males)
	output:
		directory("male_intersect.meryl")
	threads: threads,
	shell:
		"""
		meryl intersect-min threads={threads} {input.males} output {output}
		"""

#negate kmers from other sex
rule meryl_diff_F:
	input:
		F = "female_intersect.meryl",
		M = lambda wc: expand("meryl_db/{sample}.meryl", sample=males)
	output:
		directory("female-specific.meryl/"),
	shell:
		"""
		meryl difference {input.F} {input.M} output {output}
		"""
rule meryl_diff_M:
	input:
		M = "male_intersect.meryl",
		F = lambda wc: expand("meryl_db/{sample}.meryl", sample=females)
	output:
		directory("male-specific.meryl/"),
	shell:
		"""
		meryl difference {input.M} {input.F} output {output}
		"""

#remove kmer noise
rule meryl_noise_F:
	input:
		"female-specific.meryl/",
	output:
		directory("female-specific.1.meryl/"),
	shell:
		"""
		meryl greater-than 2 {input} output {output}
		"""
rule meryl_noise_M:
	input:
		"male-specific.meryl/",
	output:
		directory("male-specific.1.meryl/"),
	shell:
		"""
		meryl greater-than 2 {input} output {output}
		"""

#find sex-specific F reads
rule meryl_lookup_F1:
	input:
		r1 = "fastqs/{sample}.R1.fq.gz",
		mers = "female-specific.1.meryl/",
	output:
		"filtered_reads/females/{sample}.Fspec.R1.fastq.gz",
	shell:
		"""
		mkdir -p filtered_reads/
		mkdir -p filtered_reads/females/
		meryl-lookup -sequence {input.r1} -mers {input.mers} -include -output {output}
		"""
rule meryl_lookup_F2:
	input:
		r2 = "fastqs/{sample}.R2.fq.gz",
		mers = "female-specific.1.meryl/",
	output:
		"filtered_reads/females/{sample}.Fspec.R2.fastq.gz",
	shell:
		"""
		mkdir -p filtered_reads/
		mkdir -p filtered_reads/females/
		meryl-lookup -sequence {input.r2} -mers {input.mers} -include -output {output}
		"""

#find sex-specific M reads
rule meryl_lookup_M1:
	input:
		r1 = "fastqs/{sample}.R1.fq.gz",
		mers = "male-specific.1.meryl/",
	output:
		"filtered_reads/males/{sample}.Mspec.R1.fastq.gz",
	shell:
		"""
		mkdir -p filtered_reads/
		mkdir -p filtered_reads/males/
		meryl-lookup -sequence {input.r1} -mers {input.mers} -include -output {output}
		"""
rule meryl_lookup_M2:
	input:
		r2 = "fastqs/{sample}.R2.fq.gz",
		mers = "male-specific.1.meryl/",
	output:
		"filtered_reads/males/{sample}.Mspec.R2.fastq.gz",
	shell:
		"""
		mkdir -p filtered_reads/
		mkdir -p filtered_reads/males/
		meryl-lookup -sequence {input.r2} -mers {input.mers} -include -output {output}
		"""

#find sex-specific read pairs
rule repair_reads_F:
	input:
		r1 = "filtered_reads/females/{sample}.Fspec.R1.fastq.gz",
		r2 = "filtered_reads/females/{sample}.Fspec.R2.fastq.gz",
	output:
		r1 = "repaired_reads/females/{sample}.Fspec.R1.fq.gz",
		r2 = "repaired_reads/females/{sample}.Fspec.R2.fq.gz",
		S  = "repaired_reads/females/{sample}.Fspec.S.fastq.gz",
	resources:
		mem_gb = repair_mem,
	shell:
		"""
		mkdir -p repaired_reads/
		mkdir -p repaired_reads/females/
		repair.sh -Xmx{resources.mem_gb}g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.S} -da showspeed=t overwrite=t 2>/dev/null
		"""
rule repair_reads_M:
	input:
		r1 = "filtered_reads/males/{sample}.Mspec.R1.fastq.gz",
		r2 = "filtered_reads/males/{sample}.Mspec.R2.fastq.gz",
	output:
		r1 = "repaired_reads/males/{sample}.Mspec.R1.fq.gz",
		r2 = "repaired_reads/males/{sample}.Mspec.R2.fq.gz",
		S  = "repaired_reads/males/{sample}.Mspec.S.fastq.gz",
	resources:
		mem_gb = repair_mem,
	shell:
		"""
		mkdir -p repaired_reads/
		mkdir -p repaired_reads/males/
		repair.sh -Xmx{resources.mem_gb}g in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outs={output.S} -da showspeed=t overwrite=t 2>/dev/null
		"""

##map unique reads from in each individual
rule map_scinkd_F:
	input:
		r1 = "repaired_reads/females/{sample}.Fspec.R1.fq.gz",
		r2 = "repaired_reads/females/{sample}.Fspec.R2.fq.gz",
	output:
		bam = "mapped/females/{sample}.bam",
		bai = "mapped/females/{sample}.bam.bai",
	params:
		genome = genome,
	resources:
		mem_gb = repair_mem,
	threads: threads,
	shell:
		"""
		mkdir -p mapped/
		mkdir -p mapped/females/
		minimap2 -ax sr -t{threads} --secondary=no {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
		"""
rule map_scinkd_M:
	input:
		r1 = "repaired_reads/males/{sample}.Mspec.R1.fq.gz",
		r2 = "repaired_reads/males/{sample}.Mspec.R2.fq.gz",
	output:
		bam = "mapped/males/{sample}.bam",
		bai = "mapped/males/{sample}.bam.bai",
	params:
		genome = genome,
	threads: threads,
	resources:
		mem_gb = repair_mem,
	shell:
		"""
		mkdir -p mapped/
		mkdir -p mapped/males/
		minimap2 -ax sr -t{threads} --secondary=no {params.genome} {input.r1} {input.r2} | samtools sort --write-index -@{threads} -O bam - -o {output.bam}##idx##{output.bai} 2>/dev/null
		"""

#calculate coverage/depth for each individual
rule calc_cov_F:
	input:
		bam = "mapped/females/{sample}.bam",
	output:
		cov = "coverage/females/{sample}.regions.bed.gz",
		sum = "coverage/females/{sample}.mosdepth.summary.txt",
	params:
		prefix = lambda wc: "coverage/females/" +  wc.sample,
		window = window,
	threads: threads,
	shell:
		"""
		mkdir -p coverage/females/
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""
rule calc_cov_M:
	input:
		bam = "mapped/males/{sample}.bam",
	output:
		cov = "coverage/males/{sample}.regions.bed.gz",
		sum = "coverage/males/{sample}.mosdepth.summary.txt",
	params:
		prefix = lambda wc: "coverage/males/" +  wc.sample,
		window = window,
	threads: threads,
	shell:
		"""
		mkdir -p coverage/males/
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""

#plot uncorrected sex-specific read depth for each individual
rule plot_each_F:
	input:
		cov = "coverage/females/{sample}.regions.bed.gz",
	output:
		png = "coverage/females/{sample}.regions.manhattan.png",
		pdf = "coverage/females/{sample}.regions.manhattan.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/manhattan_plot_indiv.R {input.cov} {params.ChrNum}
		"""
rule plot_each_M:
	input:
		cov = "coverage/males/{sample}.regions.bed.gz",
	output:
		png = "coverage/males/{sample}.regions.manhattan.png",
		pdf = "coverage/males/{sample}.regions.manhattan.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/manhattan_plot_indiv.R {input.cov} {params.ChrNum}
		"""

##combine data from each individual to summarize and amplify signal
rule cat_Fspec_reads:
	input:
		r1 = expand("repaired_reads/females/{sample}.Fspec.R1.fq.gz", sample=females),
		r2 = expand("repaired_reads/females/{sample}.Fspec.R2.fq.gz", sample=females),
	output:
		cat1 = "repaired_reads/Fspec.R1.fastq.gz",
		cat2 = "repaired_reads/Fspec.R2.fastq.gz",
	shell:
		"""
		cat {input.r1} > {output.cat1}
		cat {input.r2} > {output.cat2}		
		"""
rule cat_Mspec_reads:
	input:
		r1 = expand("repaired_reads/males/{sample}.Mspec.R1.fq.gz", sample=males),
		r2 = expand("repaired_reads/males/{sample}.Mspec.R2.fq.gz", sample=males),
	output:
		cat1 = "repaired_reads/Mspec.R1.fastq.gz",
		cat2 = "repaired_reads/Mspec.R2.fastq.gz",
	shell:
		"""
		cat {input.r1} > {output.cat1}
		cat {input.r2} > {output.cat2}		
		"""

#merge M and F bam files
rule merge_bam_F:
	input:
		bam = expand("mapped/females/{sample}.bam", sample=females),
	output:
		bam = "mapped/Fspec.bam",
		bai = "mapped/Fspec.bam.bai",
	params:
		genome = genome,
	threads: threads,
	shell:
		"""
		samtools merge -@{threads} -f --write-index -o {output.bam}##idx##{output.bai} {input.bam}
		"""
rule merge_bam_M:
	input:
		bam = expand("mapped/males/{sample}.bam", sample=males),
	output:
		bam = "mapped/Mspec.bam",
		bai = "mapped/Mspec.bam.bai",
	params:
		genome = genome,
	threads: threads,
	shell:
		"""
		samtools merge -@{threads} -f --write-index -o {output.bam}##idx##{output.bai} {input.bam}
		"""

#calculate coverage/depth for combined throughput
rule cov_Fspec:
	input:
		bam = "mapped/Fspec.bam",
	output:
		cov = "coverage/Fspec_total.regions.bed.gz",
		sum = "coverage/Fspec_total.mosdepth.summary.txt",
	params:
		prefix = "coverage/Fspec_total",
		window = window,
	threads: threads,
	shell:
		"""
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""
rule cov_Mspec:
	input:
		bam = "mapped/Mspec.bam",
	output:
		cov = "coverage/Mspec_total.regions.bed.gz",
		sum = "coverage/Mspec_total.mosdepth.summary.txt",
	params:
		prefix = "coverage/Mspec_total",
		window = window,
	threads: threads,
	shell:
		"""
		MOSDEPTH_PRECISION=5 mosdepth -Q 10 -x -n --by {params.window} -t{threads} {params.prefix} {input.bam}
		"""

#plot total uncorrected read depth
rule plot_Fspec:
	input:
		cov = "coverage/Fspec_total.regions.bed.gz",
	output:
		png = "coverage/Fspec_total.regions.manhattan.png",
		pdf = "coverage/Fspec_total.regions.manhattan.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/manhattan_plot_indiv.R {input.cov} {params.ChrNum}
		"""
rule plot_Mspec:
	input:
		cov = "coverage/Mspec_total.regions.bed.gz",
	output:
		png = "coverage/Mspec_total.regions.manhattan.png",
		pdf = "coverage/Mspec_total.regions.manhattan.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/manhattan_plot_indiv.R {input.cov} {params.ChrNum}
		"""

#make dotplot
rule plot_scinkd_dotplot:
	input:
		F = "coverage/Fspec_total.mosdepth.summary.txt",
		M = "coverage/Mspec_total.mosdepth.summary.txt",
	output:
		png = "SCINKD3.dotplot.png",
		pdf = "SCINKD3.dotplot.pdf",
	params:
		ChrNum = ChrNum,
		genome = genome,
	shell:
		"""
		Rscript SCINKD3/chrom_dotplot.R SCINKD3 {input.F} {input.M} {params.ChrNum}
		"""

#plot corrected and negated read depth
rule plot_scinkd_F:
	input:
		F = "coverage/Fspec_total.regions.bed.gz",
		M = "coverage/Mspec_total.regions.bed.gz",
	output:
		png = "SCINKD3.females.png",
		pdf = "SCINKD3.females.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/SCINKD3_F-out_manhattan.R {input.M} {input.F} {params.ChrNum}
		"""
rule plot_scinkd_M:
	input:
		F = "coverage/Fspec_total.regions.bed.gz",
		M = "coverage/Mspec_total.regions.bed.gz",
	output:
		png = "SCINKD3.males.png",
		pdf = "SCINKD3.males.pdf",
	params:
		ChrNum = ChrNum,
	shell:
		"""
		Rscript SCINKD3/SCINKD3_M-out_manhattan.R {input.M} {input.F} {params.ChrNum}
		"""
