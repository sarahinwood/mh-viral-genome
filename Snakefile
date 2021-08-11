#!/usr/bin/env python3

import peppy
import pathlib2

#############
# FUNCTIONS #
#############

##needed to get BUSCO running in new folder
def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def get_peptide_dbs(wildcards):
    input_keys = ['fasta']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_extensions = pep.sample_table['sample_name']

##############
# CONTAINERS #
##############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
reapr_container = 'docker://biocontainers/reapr:v1.0.18dfsg-4-deb_cv1'
busco_container = 'docker://ezlabgva/busco:v5.2.1_cv1'
blast_container = 'docker://ncbi/blast:2.12.0'
hmmer_container = 'docker://biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
assembly_stats_container = 'docker://sangerpathogens/assembly-stats'

#########
# RULES #
#########

rule target:
    input:
        ## initial viral contigs ##
        'output/bb_stats_orig/viral_input/bb_stats.out',
        'output/assembly_stats/Mh_virus_stats.out',
        'output/hrs_self_blast/blastn.outfmt6',
        'output/hmmer/prodigalPFAM_domtblout.out',
        #'output/reapr/05.summary.report.txt',
        'output/busco/auto_prok/short_summary.specific.baculoviridae_odb10.auto_prok.txt',
        'output/interpro/Mh_interpro_table.csv',
        'output/prodigal_blast/blastp_gff.csv',
        ## PRICE output ##
        expand('output/bb_stats/{extension}_bb_stats.out', extension=all_extensions),
        expand('output/prodigal/{extension}/blastp_gff.csv', extension=all_extensions)

#####################
### PRICE contigs ###
#####################

##blastp gene predictions
rule extension_prodigal_blastp_analysis:
    input:
        prodigal_blastp_res = 'output/prodigal/{extension}/prodigal_blastp.outfmt6',
        prodigal_gff_file = 'output/prodigal/{extension}/gene_predictions.gff'
    output:
        prodigal_blastp_res_table = 'output/prodigal/{extension}/prodigal_blastp_res.csv',
        blastp_gff = 'output/prodigal/{extension}/blastp_gff.csv'
    threads:
        20
    log:
        'output/logs/{extension}/prodigal_blastp_analysis.log'
    script:
        'src/extension_prodigal_blastp_analysis.R'

rule blast_extension_prodigal_predictions:
    input:
        prodigal = 'output/prodigal/{extension}/protein_translations.faa'
    output:
        blastp_res = 'output/prodigal/{extension}/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        20
    log:
        'output/logs/{extension}/prodigal_blastp.log'
    shell:
        'blastp '
        '-query {input.prodigal} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-max_target_seqs 1 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

##re-predict viral genes on PRICE contigs using bacterial translation code
rule prodigal:
    input:
        extensions = 'joseph_output/{extension}.fasta'
    output:
        protein_translations = 'output/prodigal/{extension}/protein_translations.faa',
        nucleotide_seq = 'output/prodigal/{extension}/nucleotide_seq.fasta',
        gene_predictions = 'output/prodigal/{extension}/gene_predictions.gff'
    log:
        'output/logs/{extension}/prodigal.log'
    threads:
        1
    shell:
        'bin/Prodigal-2.6.1/prodigal '
        '-i {input.extensions} '
        '-a {output.protein_translations} '
        '-d {output.nucleotide_seq} '
        '-f gff '
        '-p meta '
        '-o {output.gene_predictions} '
        '2> {log} '


##PRICE viral contig stats
rule bb_stats_extensions:
    input:
        extensions = 'joseph_output/{extension}.fasta'
    output:
        gc = 'output/bb_stats/{extension}_gc.txt',
        stats = 'output/bb_stats/{extension}_bb_stats.out',
        gc_hist = 'output/bb_stats/{extension}_gc_hist.out'
    log:
        'output/logs/bbstats/{extension}.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.extensions} '
        'out={output.stats} '
        'gc={output.gc} '
        'gcformat=4 '
        'gchist={output.gc_hist} '
        '2> {log}'

#######################
### initial contigs ###
#######################

rule busco:
    input:
        contigs='data/Mh_DNA_virus_contigs.faa'
    output:
        'output/busco/auto_prok/short_summary.specific.baculoviridae_odb10.auto_prok.txt'
    log:
        str(pathlib2.Path(resolve_path('output/logs/'),
                            'busco.log'))
    params:
        outdir = 'auto_prok',
        wd = 'output/busco',
        contigs = lambda wildcards, input: resolve_path(input.contigs)
    singularity:
        busco_container
    threads:
        20
    shell:
        'cd {params.wd} || exit 1 ; '
        'busco '
        '--in {params.contigs} '
        '--out {params.outdir} '
        '--auto-lineage-prok '
        '--mode genome '
        '--cpu {threads} '
        '2> {log}'

##can find places contigs have been misassembled - keeps crashing
rule reapr_pipeline:
    input:
        contigs='data/Mh_DNA_virus_contigs.faa',
        smaltmap_bam = 'output/reapr/smaltmap/mapped.bam',
        perfectmap = 'output/reapr/perfect_cov.gz'
    output:
        'output/reapr/05.summary.report.txt'
    params:
        perfect_pref = 'output/reapr/',
        outdir = 'output/reapr'
    singularity:
        reapr_container
    threads:
        20
    log:
        'output/logs/reapr_pipeline.log'
    shell:
        'reapr pipeline {input.contigs} '
        '{input.smaltmap_bam} '
        '{params.outdir} '
        '{params.perfect_pref} '
        '2> {log}'

rule reapr_perfectmap:
    input:
        contigs = 'data/Mh_DNA_virus_contigs.faa',
        r1 = 'data/bbduk_trim_dna/Mh/Mh_trimr1.fq.gz',
        r2 = 'data/bbduk_trim_dna/Mh/Mh_trimr2.fq.gz'
    output:
        'output/reapr/perfect_cov.gz'
    params:
        out_pref = 'output/reapr/'
    singularity:
        reapr_container
    threads:
        20
    log:
        'output/logs/reapr_perfectmap.log'
    shell:
        'reapr perfectmap '
        '{input.contigs} '
        '{input.r1} {input.r2} '
        '426 {params.out_pref} '
        '2> {log}'

rule reapr_smaltmap:
    input:
        contigs = 'data/Mh_DNA_virus_contigs.faa',
        r1 = 'data/bbduk_trim_dna/Mh/Mh_trimr1.fq.gz',
        r2 = 'data/bbduk_trim_dna/Mh/Mh_trimr2.fq.gz'
    output:
        'output/reapr/smaltmap/mapped.bam'
    singularity:
        reapr_container
    threads:
        20
    log:
        'output/logs/reapr_smaltmap.log'
    shell:
        'reapr smaltmap '
        '{input.contigs} '
        '{input.r1} {input.r2} '
        '{output} '
        '2> {log}'

####################
## repeat regions ##
####################

rule hrs_self_blast:
    input:
        viral_contigs = 'data/Mh_DNA_virus_contigs.faa',
        db = 'output/blastdb/mh_viral_contigs.nhr'
    output:
        blast_res = 'output/hrs_self_blast/blastn.outfmt6'
    params:
        db = 'output/blastdb/mh_viral_contigs'
    threads:
        20
    log:
        'output/logs/hrs_self_blast.log'
    shell:
        'blastn '
        '-query {input.viral_contigs} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}'    

rule make_blast_db:
    input:
        'data/Mh_DNA_virus_contigs.faa'
    output:
        'output/blastdb/mh_viral_contigs.nhr'
    params:
        db_name = 'mh_viral_contigs',
        db_dir = 'output/blastdb/mh_viral_contigs'
    threads:
        10
    log:
        'output/logs/make_blast_db.log'
    shell:
        'makeblastdb '
        '-in {input} '
        '-dbtype nucl '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

#################
## basic stats ##
#################

rule assembly_stats:
    input:
        'data/Mh_DNA_virus_contigs.faa'
    output: 
        'output/assembly_stats/Mh_virus_stats.out'
    singularity:
        assembly_stats_container
    shell:
        'assembly-stats '
        '{input} > {output}'

rule bb_stats:
    input:
        viral_contigs = 'data/Mh_DNA_virus_contigs.faa'
    output:
        gc = 'output/bb_stats_orig/viral_input/gc.txt',
        stats = 'output/bb_stats_orig/viral_input/bb_stats.out',
        gc_hist = 'output/bb_stats_orig/viral_input/gc_hist.out'
    log:
        'output/logs/bbstats/viral_input.log'
    singularity:
        bbduk_container
    shell:
        'stats.sh '
        'in={input.viral_contigs} '
        'out={output.stats} '
        'gc={output.gc} '
        'gcformat=4 '
        'gchist={output.gc_hist} '
        '2> {log}'

################
## annotation ##
################

##hmmscan
rule hmmscan:
    input:
        peptides_viral_contigs = 'data/Mh_prodigal/protein_translations.faa',
        db = 'bin/db/trinotate_dbs/Pfam-A.hmm'
    output:
        'output/hmmer/prodigalPFAM_domtblout.out'
    threads:
        10
    log:
        'output/logs/hmmscan.log'
    singularity:
        hmmer_container
    shell:
        'hmmscan '
        '--cpu {threads} '
        '--domtblout {output} '
        '{input.db} '
        '{input.peptides_viral_contigs} '
        '> {log}'

##blast with lower threshold
rule prodigal_blastp_analysis:
    input:
        prodigal_blastp_res = 'output/prodigal_blast/prodigal_blastp.outfmt6',
        prodigal_gff_file = 'data/Mh_prodigal/gene_predictions.gff'
    output:
        prodigal_blastp_res_table = 'output/prodigal_blast/prodigal_blastp_res.csv',
        blastp_gff = 'output/prodigal_blast/blastp_gff.csv',
        species_no_hits = 'output/prodigal_blast/species_hit_counts.csv'
    threads:
        20
    log:
        'output/logs/prodigal_blastp_analysis.log'
    script:
        'src/prodigal_blastp_analysis.R'

rule blast_prodigal_predictions:
    input:
        prodigal = 'data/Mh_prodigal/protein_translations.faa'
    output:
        blastp_res = 'data/prodigal_blast/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blastdb/nr/nr'
    threads:
        20
    singularity:
        blast_container
    log:
        'output/logs/prodigal_blastp.log'
    shell:
        'blastp '
        '-query {input.prodigal} '
        '-db {params.blast_db} '
        '-num_threads {threads} '
        '-evalue 5e-04 '
        '-outfmt "6 std staxids salltitles" > {output.blastp_res} '
        '2> {log}'

