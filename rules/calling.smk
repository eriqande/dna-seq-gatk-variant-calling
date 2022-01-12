if "restrict-regions" in config["processing"]:

    rule compose_regions:
        input:
            config["processing"]["restrict-regions"],
        output:
            "called/{contig}.regions.bed",
        conda:
            "../envs/bedops.yaml"
        shell:
            "bedextract {wildcards.contig} {input} > {output}"


rule call_variants:
    input:
        bam=get_sample_bams,
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
        #known="resources/variation.noiupac.vcf.gz",
        #tbi="resources/variation.noiupac.vcf.gz.tbi",
        regions=(
            "called/{contig}.regions.bed"
            if config["processing"].get("restrict-regions")
            else []
        ),
    output:
        gvcf=protected("called/{sample}.{contig}.g.vcf.gz"),
    log:
        "logs/gatk/haplotypecaller/{sample}.{contig}.log",
    params:
        extra=get_call_variants_params,
    wrapper:
        "0.59.0/bio/gatk/haplotypecaller"



# this is the straight-up simple version that I use to just create
# a GVCF from the bam in mkdup, over all the regions.
rule eca_call_variants:
    input:
        bam="results/mkdup/{sample}-{unit}.bam",
        bai="results/mkdup/{sample}-{unit}.bai",
        ref="resources/genome.fasta",
        idx="resources/genome.dict",
    output:
        gvcf=protected("results/gvcf/{sample}-{unit}.g.vcf.gz"),
        idx=protected("results/gvcf/{sample}-{unit}.g.vcf.gz.tbi"),
    conda:
        "../envs/gatk4.yaml"
    log:
        stderr="results/logs/gatk/haplotypecaller/{sample}-{unit}.stderr",
        stdout="results/logs/gatk/haplotypecaller/{sample}-{unit}.stdout",
    params:
        java_opts="-Xmx4g"
    resources:
        mem_mb = 4600,
        cpus = 1
    threads: 1
    shell:
        "gatk --java-options \"{params.java_opts}\" HaplotypeCaller "
        " -R {input.ref} "
        " -I {input.bam} "
        " -O {output.gvcf} "
        " --native-pair-hmm-threads {threads} "
        " -ERC GVCF > {log.stdout} 2> {log.stderr} "


rule combine_calls:
    input:
        ref="resources/genome.fasta",
        gvcfs=expand("called/{sample}.{{contig}}.g.vcf.gz", sample=samples.index),
    output:
        gvcf="called/all.{contig}.g.vcf.gz",
    log:
        "logs/gatk/combinegvcfs.{contig}.log",
    wrapper:
        "0.59.2/bio/gatk/combinegvcfs"


rule genotype_variants:
    input:
        ref="resources/genome.fasta",
        gvcf="called/all.{contig}.g.vcf.gz",
    output:
        vcf=temp("genotyped/all.{contig}.vcf.gz"),
    params:
        extra=config["params"]["gatk"]["GenotypeGVCFs"],
    log:
        "logs/gatk/genotypegvcfs.{contig}.log",
    wrapper:
        "0.59.2/bio/gatk/genotypegvcfs"


rule merge_variants:
    input:
        vcfs=lambda w: expand("genotyped/all.{contig}.vcf.gz", contig=get_contigs()),
    output:
        vcf="genotyped/all.vcf.gz",
    log:
        "logs/picard/merge-genotyped.log",
    wrapper:
        "0.59.2/bio/picard/mergevcfs"
