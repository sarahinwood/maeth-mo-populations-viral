#!/usr/bin/env python3

import os
import peppy

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples=pep.sample_table["sample_name"]

bbduk_adapters = '/adapters.fa'
bbduk_ref = '/phix174_ill.ref.fa.gz'

##############
# CONTAINERS #
##############

tidyverse_container = 'shub://TomHarrop/singularity-containers:r_3.5.0'
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'


#########
# RULES #
#########

rule target:
    input:
        expand('output/depth_analysis/{sample}_boxplot.jpeg', sample=all_samples),
        expand('output/depth_stats/{sample}/wilcox_res.txt', sample=all_samples),
        expand('output/depth_stats_no/{sample}/wilcox_res.txt', sample=all_samples),
        'output/depth_analysis/depth_boxplot_panel.pdf',
        'output/depth_analysis/meandepth_boxplot_panel_nofacet.svg'

rule depth_stat_test_no:
    input:
        mo_gc_table = 'data/MO_contig_ids.csv',
        coverage_file = 'output/samtools_coverage/{sample}/coverage.out'
    output:
        shapiro_res = 'output/depth_stats_no/{sample}/gc_shapiro_res.txt',
        summary_stats = 'output/depth_stats_no/{sample}/summary_stats.txt',
        kruskal_res = 'output/depth_stats_no/{sample}/kruskal_res.txt',
        wilcox_res = 'output/depth_stats_no/{sample}/wilcox_res.txt'
    log:
        'output/logs/depth_stats_no/{sample}_depth_stats.log'
    script:
        'src/depth_stat_test.R'

rule depth_stat_test:
    input:
        mo_gc_table = 'data/MO_contig_ids.csv',
        coverage_file = 'output/samtools_coverage/{sample}/coverage.out'
    output:
        shapiro_res = 'output/depth_stats/{sample}/gc_shapiro_res.txt',
        summary_stats = 'output/depth_stats/{sample}/summary_stats.txt',
        kruskal_res = 'output/depth_stats/{sample}/kruskal_res.txt',
        wilcox_res = 'output/depth_stats/{sample}/wilcox_res.txt'
    log:
        'output/logs/depth_stats/{sample}_depth_stats.log'
    script:
        'src/depth_stat_test.R'

rule depth_boxplot_panel:
    input:
        st_depth_file_i44 = 'output/samtools_depth/indiv44_maethm_lincoln_depth.out',
        st_depth_file_i4 = 'output/samtools_depth/indiv4_maethm_lincoln_depth.out',
        st_depth_file_i68 = 'output/samtools_depth/indiv68_maethm_lincoln_depth.out',
        st_depth_file_i76 = 'output/samtools_depth/indiv76_maethm_lincoln_depth.out',
        mh_gc_table = 'data/MO_contig_ids.csv'
    output:
        boxplot_panel = 'output/depth_analysis/depth_boxplot_panel.pdf'
    singularity:
        tidyverse_container
    threads:
        20
    log:
        'output/logs/boxplots/depth_boxplot_panel.log'
    script:
        'src/depth_boxplot_panel.R'

rule meandepth_boxplot_panel:
    input:
        st_depth_file_i44 = 'output/samtools_coverage/indiv44_maethm_lincoln/coverage.out',
        st_depth_file_i4 = 'output/samtools_coverage/indiv4_maethm_lincoln/coverage.out',
        st_depth_file_i68 = 'output/samtools_coverage/indiv68_maethm_lincoln/coverage.out',
        st_depth_file_i76 = 'output/samtools_coverage/indiv76_maethm_lincoln/coverage.out',
        mh_gc_table = 'data/MO_contig_ids.csv'
    output:
        boxplot_panel = 'output/depth_analysis/meandepth_boxplot_panel.pdf',
        boxplot_panel_nofacet = 'output/depth_analysis/meandepth_boxplot_panel_nofacet.svg'
    singularity:
        tidyverse_container
    threads:
        20
    log:
        'output/logs/boxplots/meandepth_boxplot_panel.log'
    script:
        'src/meandepth_boxplot_panel.R'

##this is using every depth measurement - should it be mean depth instead? 2000 BUSCO contigs, 7+2 viral
rule depth_boxplot:
    input:
        samtools_depth = 'output/samtools_depth/{sample}_depth.out',
        scaffold_id_table = 'data/MO_contig_ids.csv'
    output:
        boxplot_y_zoom = 'output/depth_analysis/{sample}_boxplot_y_zoom.pdf',
        boxplot = 'output/depth_analysis/{sample}_boxplot.jpeg'
    singularity:
        tidyverse_container
    threads:
        20
    log:
        'output/logs/boxplots/{sample}_depth_boxplot.log'
    script:
        'src/depth_boxplot.R'

####################
## depth/coverage ##
####################

rule samtools_coverage:
    input:
        bam = 'output/samtools/{sample}/sorted.bam'
    output:
        coverage_out = 'output/samtools_coverage/{sample}/coverage.out'
    log:
        'output/logs/samtools_coverage_{sample}.log'
    threads:
        20
    shell:
        'samtools coverage '
        '{input.bam} '
        '-o {output.coverage_out} '
        '2> {log}'

rule samtools_depth:
    input:
        sorted_bam = 'output/samtools/{sample}/sorted.bam'
    output:
        depth_out = 'output/samtools_depth/{sample}_depth.out'
    log:
        'output/logs/samtools_depth/{sample}.log'
    threads:
        20
    shell:
        'samtools depth '
        '{input.sorted_bam} '
        '-a '
        '> {output.depth_out} '
        '2> {log}'

rule samtools_sort:
    input:
        sam = 'output/bwa/{sample}/bwa_mem.sam'
    output:
        sorted_bam = 'output/samtools/{sample}/sorted.bam'
    log:
        'output/logs/samtools_sort_{sample}.log'
    threads:
        20
    shell:
        'samtools sort '
        '{input.sam} '
        '-o {output.sorted_bam} '
        '2> {log}'

##bwa-mem to map for read depth of scaffolds
rule bwa_mem:
    input:
        index = 'output/bwa/index.bwt',
        trim_fastq = 'output/bbduk_trim_dna/{sample}/{sample}_trim.fq.gz'
    output:
        sam = 'output/bwa/{sample}/bwa_mem.sam'
    params:
        index_dir = 'output/bwa/index'
    threads:
        50
    log:
        'output/logs/bwa_mem_{sample}.log'
    shell:
        'bwa mem '
        '-t {threads} '
        '{params.index_dir} '
        '{input.trim_fastq} '
        '-p '
        '> {output.sam} '
        '2> {log}'

rule bwa_index:
    input:
        genome = 'data/MO_final_genome.fa'
    output:
        index = 'output/bwa/index.bwt'
    params:
        outdir = 'output/bwa/index'
    threads:
        20
    log:
        'output/logs/bwa_index.log'
    shell:
        'bwa index '
        '{input.genome} '
        '-p {params.outdir} '
        '2> {log} '

##trim and decontaminate DNA reads to map onto genome
rule bbduk_trim_dna:
    input:
        fastq = 'output/bbduk_trim_dna/{sample}/{sample}_filr.fq.gz'
    output:
        trim_fastq = 'output/bbduk_trim_dna/{sample}/{sample}_trim.fq.gz',
        t_stats = 'output/bbduk_trim_dna/{sample}/trim-stats.txt'
    log:
        trim = 'output/logs/bbduk_trim_dna/{sample}_trim.log'
    params:
        adapters = bbduk_adapters
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.fastq} '
        'int=t '
        'out={output.trim_fastq} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo '
        'forcetrimmod=5 '
        'stats={output.t_stats} '
        '2> {log.trim} '

rule bbduk_filter_dna:  
    input:
        fastq = 'data/mo_fastq/{sample}.fastq'
    output:
        fastq = 'output/bbduk_trim_dna/{sample}/{sample}_filr.fq.gz',
        f_stats = 'output/bbduk_trim_dna/{sample}/filter-stats.txt'
    log:
        filter = 'output/logs/bbduk_trim_dna/{sample}_filter.log'
    params:
        ref = bbduk_ref
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'threads={threads} '
        'in={input.fastq} '
        'int=t '
        'out={output.fastq} '
        'ref={params.ref} '
        'hdist=1 '
        'stats={output.f_stats} '       
        '2> {log.filter} ' 