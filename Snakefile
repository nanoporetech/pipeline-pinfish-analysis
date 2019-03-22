
import os
from os import path

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)

include: "snakelib/utils.snake"

rule dump_versions:
    output:
        ver = "versions.txt"
    conda: "env.yml"
    shell:"""
    command -v conda > /dev/null && conda list > {output.ver}
    echo -e "pinfish/spliced_bam2gff\t" `{SNAKEDIR}/pinfish/spliced_bam2gff/spliced_bam2gff -V` >> {output.ver}
    echo -e "pinfish/cluster_gff\t" `{SNAKEDIR}/pinfish/cluster_gff/cluster_gff -V` >> {output.ver}
    echo -e "pinfish/collapse_partials" `{SNAKEDIR}/pinfish/collapse_partials/collapse_partials -V` >> {output.ver}
    echo -e "pinfish/polish_clusters" `{SNAKEDIR}/pinfish/polish_clusters/polish_clusters -V` >> {output.ver}
    """

rule build_minimap_index: ## build minimap2 index
    input:
        genome = config["genome_fasta"]
    output:
        index = "index/genome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """

rule map_reads: ## map reads using minimap2
    input:
       index = rules.build_minimap_index.output.index,
       fastq = config["reads_fastq"]
    output:
       bam = "alignments/reads_aln_sorted.bam"
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

rule convert_bam: ## convert BAM to GFF
    input:
        bam = rules.map_reads.output.bam
    output:
        raw_gff = "results/raw_transcripts.gff"
    conda: "env.yml"
    params:
        opts = config["spliced_bam2gff_opts"]
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/spliced_bam2gff/spliced_bam2gff {params.opts} -t {threads} -M {input.bam} > {output.raw_gff}
    """

rule cluster_gff: ## cluster transcripts in GFF
    input:
        raw_gff = rules.convert_bam.output.raw_gff
    output:
        cls_gff = "results/clustered_transcripts.gff",
        cls_tab = "results/cluster_memberships.tsv",
    params:
        c = config["minimum_cluster_size"],
        d = config["exon_boundary_tolerance"],
        e = config["terminal_exon_boundary_tolerance"],
        min_iso_frac = config["minimum_isoform_percent"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/cluster_gff/cluster_gff -p {params.min_iso_frac} -t {threads} -c {params.c} -d {params.d} -e {params.e} -a {output.cls_tab} {input.raw_gff} > {output.cls_gff}
    """

rule collapse_clustered: ## collapse clustered read artifacts
    input:
        cls_gff = rules.cluster_gff.output.cls_gff
    output:
        cls_gff_col = "results/clustered_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/pinfish/collapse_partials/collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.cls_gff} > {output.cls_gff_col}
    """

rule polish_clusters: ## polish read clusters
    input:
        cls_gff = rules.cluster_gff.output.cls_gff,
        cls_tab = rules.cluster_gff.output.cls_tab,
        bam = rules.map_reads.output.bam,
    output:
        pol_trs = "results/polished_transcripts.fas",
    params:
        c = config["minimum_cluster_size"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/polish_clusters/polish_clusters -t {threads} -a {input.cls_tab} -c {params.c} -o {output.pol_trs} {input.bam}
    """

rule map_polished: ## map polished transcripts to genome
    input:
       index = rules.build_minimap_index.output.index,
       fasta = rules.polish_clusters.output.pol_trs,
    output:
       pol_bam = "alignments/polished_reads_aln_sorted.bam"
    params:
        extra = config["minimap2_opts_polished"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} {params.extra} -ax splice {input.index} {input.fasta}\
    | samtools view -Sb -F 2304 | samtools sort -@ {threads} - -o {output.pol_bam};
    samtools index {output.pol_bam}
    """

rule convert_polished: ## convert BAM of polished transcripts to GFF
    input:
        bam = rules.map_polished.output.pol_bam
    output:
        pol_gff = "results/polished_transcripts.gff"
    params:
        extra = config["spliced_bam2gff_opts_pol"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/pinfish/spliced_bam2gff/spliced_bam2gff {params.extra} -t {threads} -M {input.bam} > {output.pol_gff}
    """

rule collapse_polished: ## collapse polished read artifacts
    input:
        pol_gff = rules.convert_polished.output.pol_gff
    output:
        pol_gff_col = "results/polished_transcripts_collapsed.gff"
    params:
        d = config["collapse_internal_tol"],
        e = config["collapse_three_tol"],
        f = config["collapse_five_tol"],
    conda: "env.yml"
    shell:"""
    {SNAKEDIR}/pinfish/collapse_partials/collapse_partials -d {params.d} -e {params.e} -f {params.f} {input.pol_gff} > {output.pol_gff_col}
    """

rule gen_corr_trs: ## Generate corrected transcriptome.
    input:
        genome = config["genome_fasta"],
        gff = rules.collapse_polished.output.pol_gff_col,
    output:
        fasta = "results/corrected_transcriptome_polished_collapsed.fas"
    conda: "env.yml"
    shell:"""
    gffread -g {input.genome} -w {output.fasta} {input.gff}
    """

rule all: ## run the whole pipeline
    input:
        ver = rules.dump_versions.output.ver,
        index = rules.build_minimap_index.output.index,
        aligned_reads = rules.map_reads.output.bam,
        raw_gff = rules.convert_bam.output.raw_gff,
        cls_gff = rules.cluster_gff.output.cls_gff,
        cls_gff_col = rules.collapse_clustered.output.cls_gff_col,
        pol_trs = rules.polish_clusters.output.pol_trs,
        pol_bam = rules.map_polished.output.pol_bam,
        pol_gff = rules.convert_polished.output.pol_gff,
        pol_gff_col = rules.collapse_polished.output.pol_gff_col,
        corr_trs = rules.gen_corr_trs.output.fasta,

