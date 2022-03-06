configfile: "metadata/config.yaml"

import pandas as pd
import os

units_table = pd.read_table(config["samples_file"], sep="\t").set_index("Run", drop=False)
# samples= list(units_table.source_name.unique())

groups=pd.read_csv(config["comparisons"], sep ="\t")
# control_groups=list(groups.Control_Group)
# treatment_groups=list(groups.Treatment_Group)

rule all:
	input:
		expand('files/deltapsi_files/{comparison}.deltapsi.tsv',comparison=groups.Comparison),
		expand("files/bam_files/{sample}_Aligned.sortedByCoord.out.bam", sample=units_table.Sample_Name),
		expand("results/salmon/{unit}_quant/quant.sf",unit=units_table.Run)


rule SRR_download:
	params:
		srr=lambda wildcards: units_table.Run[units_table.Run == wildcards.unit],
		path="files/fastq_files/"
	output:
		"files/fastq_files/{unit}/{unit}.sra"
	conda:
		"envs/SraTools.yaml"
	resources:
		cpu=1,
		mem=lambda wildcards, attempt: attempt * 4,
		time=360
	shell:
		"prefetch -O {params.path} {wildcards.unit}"

rule SRR_convert:
	input:
		srr="files/fastq_files/{unit}/{unit}.sra"
	params:
		path="files/fastq_files/"
	output:
		expand("files/fastq_files/{{unit}}_{num}.fastq.gz", num=[1,2])
	conda:
		"envs/SraTools.yaml"
	resources:
		cpu=1,
		mem=lambda wildcards, attempt: attempt * 4,
		time=360
	shell:
		"fastq-dump --gzip --split-3 -O {params.path} {input.srr}"



rule STAR_align:
	input:
		genome=config["STAR_genome"],
		files=lambda wildcards: expand("files/fastq_files/{unit}_{num}.fastq.gz", unit=units_table.Run[units_table.Sample_Name == wildcards.sample], num=[1,2])
	params:
		path="files/bam_files/{sample}_"
	resources:
		cpu=10,
		mem=lambda wildcards, attempt: attempt * 120 
	output:
		"files/bam_files/{sample}_Aligned.sortedByCoord.out.bam"
	shell:
		"STAR --twopassMode Basic --genomeDir {input.genome} --outTmpKeep None "
		"--readFilesIn {input.files} --readFilesCommand zcat "
		"--runThreadN {resources.cpu} --outSAMtype BAM SortedByCoordinate "
		"--outFileNamePrefix {params.path} --alignSJoverhangMin 8 "
		"--limitBAMsortRAM {resources.mem}000000000 --outSAMattributes All "
		"--quantMode GeneCounts"

# rule samtools_merge_bam:
# 	input:
# 		lambda wildcards: expand('files/bam_files/{unit}_Aligned.sortedByCoord.out.bam',
# 			unit= units_table.Run[units_table.source_name == wildcards.sample])
# 	output:
# 		bam = 'files/bam_files/{sample}_merged.bam'
# 	resources:
# 		cpu=2,
# 		mem=lambda wildcards, attempt: attempt * 8
# 	shell:  
# 		"""
# 		PS1=dummy
# 		. $(conda info --base)/etc/profile.d/conda.sh
# 		conda activate samtools
# 		samtools merge {output.bam} {input}
# 		"""
		

rule samtools_index:
	input:
		'files/bam_files/{sample}_Aligned.sortedByCoord.out.bam'
	output:
		'files/bam_files/{sample}_Aligned.sortedByCoord.out.bam.bai'
	conda:
		"envs/SamTools_env.yaml"
	resources:
		cpu=10,
		mem=lambda wildcards, attempt: attempt * 10
	shell:
		"samtools index -@ {resources.cpu} {input} {output}"


rule config_majiq:
	input:
		bam=expand('files/bam_files/{sample}_Aligned.sortedByCoord.out.bam',sample=units_table.Sample_Name),
		bai=expand('files/bam_files/{sample}_Aligned.sortedByCoord.out.bam.bai',sample=units_table.Sample_Name)
	output:
		'metadata/majiq_config_file.txt'
	resources:
		cpu=1,
		mem=lambda wildcards, attempt: attempt * 1
	params:
		samples=list(units_table.Sample_Name),
		genome="hg38"
	script:
		'scripts/config_creator_majiq.py' 

rule majiq_build:
	input:
		bam=expand('files/bam_files/{sample}_Aligned.sortedByCoord.out.bam',sample=units_table.Sample_Name),
		bai=expand('files/bam_files/{sample}_Aligned.sortedByCoord.out.bam.bai',sample=units_table.Sample_Name),
		genome=config["genome"],
		majiq_config=config["build_config"]
	params:
		path='files/majiq_files/'
	conda:
		"envs/majiq_env.yaml"
	resources:
		cpu=8,
		mem=lambda wildcards, attempt: attempt * 8
	output:
		majiq_files=expand('files/majiq_files/{sample}_Aligned.sortedByCoord.out.majiq',sample=units_table.Sample_Name),
		sj_files=expand('files/majiq_files/{sample}_Aligned.sortedByCoord.out.sj',sample=units_table.Sample_Name),
		splicegraph='files/majiq_files/splicegraph.sql'
#		log='files/majiq_files/majiq.log'
	shell:
		"majiq build --conf {input.majiq_config} --nproc {resources.cpu} "
		"--disable-ir --simplify -o {params.path} {input.genome}" 
	
# #rule to create new column with underscore between names and loop
# #this will in sense merge the two columns

rule delta_psi:
	input:
		control_majiq_files=lambda wildcards: expand(
			'files/majiq_files/{group_control}{rep}_Aligned.sortedByCoord.out.majiq',
			group_control=groups.Control_Group[groups.Comparison == wildcards.comparison],rep=["_rep1","_rep2","_rep3"]),
		treatment_majiq_files=lambda wildcards: expand(
			'files/majiq_files/{group_treatment}{rep}_Aligned.sortedByCoord.out.majiq',
			group_treatment=groups.Treatment_Group[groups.Comparison == wildcards.comparison],rep=["_rep1","_rep2","_rep3"])
	conda:
		"envs/majiq_env.yaml"
	resources:
		cpu=4,
		mem=lambda wildcards, attempt: attempt * 64
	params:
		path="files/deltapsi_files",
		control=lambda wildcards: expand(
			'{group_control}',group_control=groups.Control_Group[groups.Comparison == wildcards.comparison]),
		treated=lambda wildcards: expand(
			'{group_treatment}',group_treatment=groups.Treatment_Group[groups.Comparison == wildcards.comparison])
	output:
		'files/deltapsi_files/{comparison}.deltapsi.voila'
	shell:
		"majiq deltapsi --output-type voila -grp1 {input.control_majiq_files} -grp2 {input.treatment_majiq_files} "
		"-n {params.control} {params.treated} -j {resources.cpu} -o {params.path}"
		
rule voila_tsv:
	input:
		voila='files/deltapsi_files/{comparison}.deltapsi.voila',
		sql='files/majiq_files/splicegraph.sql'
	resources:
		cpu=4,
		mem=lambda wildcards, attempt: attempt * 64
	conda:
		"envs/majiq_env.yaml"
	output:
		'files/deltapsi_files/{comparison}.deltapsi.tsv'
	shell:
		"voila tsv -f {output} --show-all -j {resources.cpu} {input.voila} {input.sql}"

#rule psi:
#	input: 
#		majiq_files=expand('files/majiq_files/{sample}_Aligned.sortedByCoord.out.majiq',
#			sample=samples)
#	resources:
#		cpu=4,
#		mem="4G"
#	conda:
#		"envs/majiq_env.yaml"
#	params:
#		path="files/psi_files"
#	output:
#		tsv_files='files/psi_files/{sample}.psi.tsv',
#		voila_files='files/psi_files/{sample}.psi.voila'
#	shell:
#		"majiq psi --output-type all -j {resources.cpu} {input.majiq_files} -o {params.path}"



rule salmon_quant:
        input:
                "files/fastq_files/{unit}_1.fastq.gz", "files/fastq_files/{unit}_2.fastq.gz"
        resources:
                cpu=8,
                time="96:00:00",
                mem=lambda wildcards, attempt: attempt * 25
        params:
                genome=config["salmon_genome"]
        output:
                "results/salmon/{unit}_quant/quant.sf"
        shell:
                "/home/torresdizm/bin/salmon quant "
                "-i {params.genome} "
                "-l A "
                "-1 {input[0]} -2 {input[1]} "
                "-p 8 --validateMapping "
                "-o ./results/salmon/{wildcards.unit}_quant"
