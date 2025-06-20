configfile: "config.yaml"

import glob


def is_paired(wildcards):
    dir = checkpoints.convert_sra.get(sample=wildcards.sample).output[0]
    return bool(
        glob.glob(f"{dir}/{wildcards.sample}_1.fastq") and
        glob.glob(f"{dir}/{wildcards.sample}_2.fastq")
    )

def get_accessions(config):
    accession_file = config["accession_file"]
    with open(accession_file) as f:
        return [line.strip() for line in f if line.strip()]

accessions = get_accessions(config)
references = config["reference"]
eve_files = config["eves"]
index_files = ["bowtie2_index/index.{}.bt2".format(ext) for ext in [1,2,3,4,"rev.1","rev.2"]]

if config["reference"]:
    ruleorder: map_reads_bowtie2_paired > skip_map_paired
    ruleorder: map_reads_bowtie2_single > skip_map_single

    if config["eves"]:
        ruleorder: index_masked_ref > index_host_bowtie2
    else:
        ruleorder: index_host_bowtie2 > index_masked_ref
else:
    ruleorder: skip_map_paired > map_reads_bowtie2_paired
    ruleorder: skip_map_single > map_reads_bowtie2_single


rule all:
    input:
        expand("diamond_results/{sample}.tsv", sample=accessions)


checkpoint convert_sra:
    input: 
        "SRA_folder/{sample}/{sample}.sra"
    output:
        directory("FASTQ_folder/{sample}")
    params:
        "SRA_folder/{sample}"
    shell:
        """
        (
          cd {params}
          fasterq-dump --split-3 {wildcards.sample}.sra
        )
        mkdir -p {output}
        mv {params}/{wildcards.sample}*.fastq {output}
        """

rule prefetch_accession:
    output:
        sra=temp("SRA_folder/{sample}/{sample}.sra")
    shell:
        """
        prefetch {wildcards.sample} -O SRA_folder/
        """

rule map_eves_to_ref:
    params:
        eve_files=config["eves"],
        references=config["reference"]
    output:
        paf_file=temp("paf/eve.paf")
    threads: 4
    shell:
        """
        minimap2 -x asm5 -t {threads} {params.references} {params.eve_files} > {output.paf_file}
        """

rule convert_paf_to_bed:
    input:
        paf_file="paf/eve.paf"
    output:
        bed_file=temp("bed/eve.bed")
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{print $6, $8, $9, $1, $12, $5}}' {input.paf_file} > {output.bed_file}
        """

rule mask_ref:
    input:
        bed_file="bed/eve.bed"
    params:
    	references=config["reference"]
    output:
        masked_ref="masked_reference/masked_ref.fna"
    shell:
        """
        bedtools maskfasta -fi {params.references} -bed {input.bed_file} -fo {output.masked_ref} -fullHeader
        """
rule index_masked_ref:
    input: 
       ref="masked_reference/masked_ref.fna"
    output:
        index=index_files
    log:
        "logs/index_masked_bowtie2.log"
    shell:
        """
        mkdir -p logs
        echo "Building index with references: {input.ref}" > {log}
        bowtie2-build {input.ref} bowtie2_index/index >> {log} 2>&1
        """
        
rule index_host_bowtie2:
    params:
        references=config["reference"]
    output:
        index=index_files
    log:
        "logs/index_host_bowtie2.log"
    shell:
        """
        mkdir -p logs
        echo "Building index with references: {params.references}" > {log}
        bowtie2-build {params.references} bowtie2_index/index >> {log} 2>&1
        """

rule map_reads_bowtie2_paired:
    input:
        trimmed_forward="trimmed_paired/{sample}_1.fq",
        trimmed_reverse="trimmed_paired/{sample}_2.fq",
        index_files=expand("bowtie2_index/index.{ext}.bt2", ext=[1,2,3,4,"rev.1","rev.2"])
    output:
        sam=temp("bowtie2_aligned_paired/{sample}.sam"),
        fq1=temp("bowtie2_aligned_paired/{sample}_unconc.1.fq"),
        fq2=temp("bowtie2_aligned_paired/{sample}_unconc.2.fq")
    params:
        unc="bowtie2_aligned_paired/{sample}_unconc.fq",
        index_prefix="bowtie2_index/index"
    threads: 30
    shell:
        """
        mkdir -p bowtie2_aligned_paired
        bowtie2 --very-sensitive \
                --un-conc {params.unc} \
                -x {params.index_prefix} \
                -1 {input.trimmed_forward} \
                -2 {input.trimmed_reverse} \
                -p {threads} \
                -S {output.sam}
        """

rule skip_map_paired:
    input:
        trimmed_forward="trimmed_paired/{sample}_1.fq",
        trimmed_reverse="trimmed_paired/{sample}_2.fq"
    output:
        fq1=temp("bowtie2_aligned_paired/{sample}_unconc_1.fq"),
        fq2=temp("bowtie2_aligned_paired/{sample}_unconc_2.fq"),
        sam=temp("bowtie2_aligned_paired/{sample}.sam")
    shell:
        """
        cp {input.trimmed_forward} {output.fq1}
        cp {input.trimmed_reverse} {output.fq2}
        touch {output.sam}
        """

rule skip_map_single:
    input:
        trimmed_single="trimmed_single/{sample}.fq"
    output:
        single_skip=temp("bowtie2_aligned_single/{sample}_unconc.fq"),
        sam=temp("bowtie2_aligned_single/{sample}.sam")
    shell:
        """
        cp {input.trimmed_single} {output.single_skip}
        touch {output.sam}
        """

rule trim_reads_fastp_paired:
    input:
        forward="FASTQ_folder/{sample}/{sample}_1.fastq",
        reverse_reads="FASTQ_folder/{sample}/{sample}_2.fastq"
    output:
        trimmed_forward=temp("trimmed_paired/{sample}_1.fq"),
        trimmed_reverse=temp("trimmed_paired/{sample}_2.fq"),
        report_html=temp("trimmed_paired/{sample}_report.html"),
        report_json=temp("trimmed_paired/{sample}_report.json")
    params:
        quality=20,
        min_length=25,
        adapter_trim="--detect_adapter_for_pe"
    threads: 4
    shell:
        """
        fastp -i {input.forward} -I {input.reverse_reads} \
              -o {output.trimmed_forward} -O {output.trimmed_reverse} \
              -q {params.quality} -l {params.min_length} --thread {threads} \
              {params.adapter_trim} \
              --html {output.report_html} --json {output.report_json}
        """

rule map_reads_bowtie2_single:
    input:
        trimmed_single="trimmed_single/{sample}.fq",
        index=expand("bowtie2_index/index.{ext}.bt2", ext=[1,2,3,4,"rev.1","rev.2"])
    output:
        sam=temp("bowtie2_aligned_single/{sample}.sam"),
        fq_s=temp("bowtie2_aligned_single/{sample}_unconc.fq")
    params:
        unc="bowtie2_aligned_single/{sample}_unconc.fq"
    threads: 30
    shell:
        """
        bowtie2 --very-sensitive \
            --un {params.unc} \
            -x bowtie2_index/index \
            -U {input.trimmed_single} \
            -p {threads} \
            -S {output.sam}
        """

rule trim_reads_fastp_single:
    input:
        single="FASTQ_folder/{sample}/{sample}.fastq"
    output:
        trimmed_single=temp("trimmed_single/{sample}.fq"),
        report_html=temp("trimmed_single/{sample}_report.html"),
        report_json=temp("trimmed_single/{sample}_report.json")
    params:
        quality=20,
        min_length=25,
        adapter_trim="--detect_adapter_for_pe"
    threads: 4
    shell:
        """
        fastp -i {input.single} \
              -o {output.trimmed_single} \
              -q {params.quality} -l {params.min_length} --thread {threads} \
              {params.adapter_trim} \
              --html {output.report_html} --json {output.report_json}
        """

rule megahit_assembly_single:
    input:
        single_reads="bowtie2_aligned_single/{sample}_unconc.fq"
    output:
        "megahit_results_single/{sample}/final.contigs.fa"
    log:
        "logs/megahit_single/{sample}.log"
    benchmark:
        "benchmarks/megahit_single/{sample}.benchmark.txt"
    threads: 30
    shell:
        """
        rm -rf megahit_results_single/{wildcards.sample}
        megahit -r {input.single_reads} \
                -o megahit_results_single/{wildcards.sample} \
                -t {threads} --verbose > {log} 2>&1
        """

rule megahit_assembly_paired:
    input:
        forward="bowtie2_aligned_paired/{sample}_unconc.1.fq",
        reverse_reads="bowtie2_aligned_paired/{sample}_unconc.2.fq"
    output:
        "megahit_results_paired/{sample}/final.contigs.fa"
    log:
        "logs/megahit_paired/{sample}.log"
    threads: 1
    shell:
        """
        rm -rf megahit_results_paired/{wildcards.sample}
        megahit -1 {input.forward} -2 {input.reverse_reads} \
                -o megahit_results_paired/{wildcards.sample} \
                -t {threads} --verbose > {log} 2>&1
        """

rule diamond:
    input:
        branch(
            is_paired,
            then="megahit_results_paired/{sample}/final.contigs.fa",
            otherwise="megahit_results_single/{sample}/final.contigs.fa"
        )
    output:
        "diamond_results/{sample}.tsv"
    threads: 16
    shell:
        """
        diamond blastx --query {input} --db {config[diamond_db]} \
        --out {output} --outfmt 6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore \
        --threads {threads}
        """

