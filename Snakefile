#!/usr/bin/env python3

import peppy

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_fastas = pep.sample_table['filename']
all_extensions = pep.sample_table['folder']

##############
# CONTAINERS #
##############

bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
reapr_container = 'docker://biocontainers/reapr:v1.0.18dfsg-4-deb_cv1'
busco_container = 'docker://ezlabgva/busco:v5.2.1_cv1'
hmmer_container = 'docker://biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
assembly_stats_container = 'docker://sangerpathogens/assembly-stats'
prodigal_container = 'docker://biocontainers/prodigal:v1-2.6.3-4-deb_cv1'
blast_container = 'docker://ncbi/blast:2.12.0'

#########
# RULES #
#########

rule target:
    input:
        ## PRICE output ##
        expand('output/bb_stats/{fasta}_bb_stats.out', fasta=all_fastas),
        expand('output/prodigal/{fasta}/blastp_gff.csv', fasta=all_fastas),
        'output/blast_c4e30_initial/blastn.outfmt6',
        ##hrs blast
        'output/hrs_self_blast/blastn.outfmt6',
        ##Spades
        'output/prodigal/spades/prodigal_blastp.outfmt6'

####################
## repeat regions ##
####################

rule hrs_self_blast:
    input:
        viral_contigs = 'data/price2/extended4.cycle30.fasta',
        db = 'output/blastdb/price2_e4c30.nhr'
    output:
        blast_res = 'output/hrs_self_blast/blastn.outfmt6'
    params:
        db = 'output/blastdb/price2_e4c30'
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

rule make_blast_db_e4c30:
    input:
        'data/price2/extended4.cycle30.fasta'
    output:
        'output/blastdb/price2_e4c30.nhr'
    params:
        db_name = 'price2_e4c30',
        db_dir = 'output/blastdb/price2_e4c30'
    threads:
        10
    log:
        'output/logs/make_blast_db_e4c30.log'
    shell:
        'makeblastdb '
        '-in {input} '
        '-dbtype nucl '
        '-title {params.db_name} '
        '-out {params.db_dir} '
        '-parse_seqids '
        '2> {log}'

## blast e4c30 genes against original prodigal genes ##
rule blast_c4e30_initial:
    input:
        c4e30_genes = 'output/prodigal/price2/extended4.cycle30/nucleotide_seq.fasta',
        db = 'output/blastdb/initial_mh_prodigal.nhr'
    output:
        blast_res = 'output/blast_c4e30_initial/blastn.outfmt6'
    params:
        db = 'output/blastdb/initial_mh_prodigal'
    threads:
        20
    log:
        'output/logs/blast_c4e30_initial.log'
    shell:
        'blastn '
        '-query {input.c4e30_genes} '
        '-db {params.db} '
        '-num_threads {threads} '
        '-evalue 1e-05 '
        '-outfmt "6 std salltitles" > {output.blast_res} '
        '2>{log}' 

rule make_blast_db:
    input:
        'data/Mh_prodigal/nucleotide_seq.fasta'
    output:
        'output/blastdb/initial_mh_prodigal.nhr'
    params:
        db_name = 'initial_mh_prodigal',
        db_dir = 'output/blastdb/initial_mh_prodigal'
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

#####################
### PRICE contigs ###
#####################

##blastp gene predictions
rule extension_prodigal_blastp_analysis:
    input:
        prodigal_blastp_res = 'output/prodigal/{extension}/{fasta}/prodigal_blastp.outfmt6',
        prodigal_gff_file = 'output/prodigal/{extension}/{fasta}/gene_predictions.gff'
    output:
        prodigal_blastp_res_table = 'output/prodigal/{extension}/{fasta}/prodigal_blastp_res.csv',
        blastp_gff = 'output/prodigal/{extension}/{fasta}/blastp_gff.csv'
    threads:
        20
    log:
        'output/logs/{extension}/{fasta}/prodigal_blastp_analysis.log'
    script:
        'src/extension_prodigal_blastp_analysis.R'

rule blast_prodigal_price:
    input:
        prodigal = 'output/prodigal/{extension}/{fasta}/protein_translations.faa'
    output:
        blastp_res = 'output/prodigal/{extension}/{fasta}/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blast_db/nr/nr'
    threads:
        20
    log:
        'output/logs/{extension}/{fasta}/prodigal_blastp.log'
    singularity:
        blast_container
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
rule prodigal_price:
    input:
        extensions = 'data/{extension}/{fasta}.fasta'
    output:
        protein_translations = 'output/prodigal/{extension}/{fasta}/protein_translations.faa',
        nucleotide_seq = 'output/prodigal/{extension}/{fasta}/nucleotide_seq.fasta',
        gene_predictions = 'output/prodigal/{extension}/{fasta}/gene_predictions.gff'
    log:
        'output/logs/{extension}/{fasta}/prodigal.log'
    threads:
        1
    singularity:
        prodigal_container
    shell:
        'prodigal '
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
        extensions = 'data/{extension}/{fasta}.fasta'
    output:
        gc = 'output/bb_stats/{extension}/{fasta}_gc.txt',
        stats = 'output/bb_stats/{extension}/{fasta}_bb_stats.out',
        gc_hist = 'output/bb_stats/{extension}/{fasta}_gc_hist.out'
    log:
        'output/logs/{extension}/{fasta}_bbstats.log'
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

##############
### spades ###
##############

rule blast_prodigal_spades:
    input:
        prodigal = 'output/prodigal/spades/protein_translations.faa'
    output:
        blastp_res = 'output/prodigal/spades/prodigal_blastp.outfmt6'
    params:
        blast_db = 'bin/db/blast_db/nr/nr'
    threads:
        20
    log:
        'output/logs/spades/prodigal_blastp.log'
    singularity:
        blast_container
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
rule prodigal_spades:
    input:
        extensions = 'data/spades_assembly/scaffolds.fasta'
    output:
        protein_translations = 'output/prodigal/spades/protein_translations.faa',
        nucleotide_seq = 'output/prodigal/spades/nucleotide_seq.fasta',
        gene_predictions = 'output/prodigal/spades/gene_predictions.gff'
    log:
        'output/logs/spades/prodigal.log'
    threads:
        1
    singularity:
        prodigal_container
    shell:
        'prodigal '
        '-i {input.extensions} '
        '-a {output.protein_translations} '
        '-d {output.nucleotide_seq} '
        '-f gff '
        '-p meta '
        '-o {output.gene_predictions} '
        '2> {log} '