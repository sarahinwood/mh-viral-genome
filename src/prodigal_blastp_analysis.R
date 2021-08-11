#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(dplyr)
library(rtracklayer)

###########
# GLOBALS #
###########

prodigal_blastp_res <- snakemake@input[["prodigal_blastp_res"]]
prodigal_gff_file <- snakemake@input[["prodigal_gff_file"]]

########
# MAIN #
########

prodigal_gff <- readGFF(prodigal_gff_file)
prodigal_blastp <- fread(prodigal_blastp_res)

setnames(prodigal_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("prodigal_nt_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##remove hypothetical protein from bacterium hits
prodigal_blastp <- subset(prodigal_blastp, taxid!=1869227)
##remove pieris annotations
prodigal_blastp <- subset(prodigal_blastp, taxid!=345717)

##sum of hits for each gene
no_hits <- count(prodigal_blastp, prodigal_nt_id)

##species hits table
prodigal_blastp_species <- prodigal_blastp
prodigal_blastp_species$species <- tstrsplit(prodigal_blastp_species$annotation, "[", fixed=TRUE, keep=c(2))
prodigal_blastp_species$species <- tstrsplit(prodigal_blastp_species$species, "]", fixed=TRUE, keep=c(1))
##counts table
species_no_hits <- count(prodigal_blastp_species, species, prodigal_nt_id)
species_no_hits_virus <- dplyr::filter(species_no_hits, grepl("virus", species))
species_gene_no_hits <- count(species_no_hits_virus, species)
fwrite(species_gene_no_hits, snakemake@output[["species_no_hits"]])

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(prodigal_blastp, prodigal_nt_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- prodigal_blastp[,.SD[which.min(evalue)], by=prodigal_nt_id]

no_hits_best_hit <- merge(no_hits, min_evalues)
no_hits_annots <- no_hits_best_hit[,c(1,2,15,12,4)]
fwrite(no_hits_annots, "output/prodigal_blast/no_hits_by_best_annot.csv")

IAP_hits <- subset(prodigal_blastp, prodigal_blastp$prodigal_nt_id=="scaffold_70_11")
bacterium_hits <- subset(prodigal_blastp, prodigal_blastp$prodigal_nt_id=="scaffold_1039_2")
bacterium2_hits <- subset(prodigal_blastp, prodigal_blastp$prodigal_nt_id=="scaffold_1195_1")

##merge blastp res with prodigal
gene_coords <- subset(prodigal_gff, select=c(start, end, ID, conf, partial, strand))
##scaffold id to gene id
scaffold_to_geneid <- data.table(prodigal_gff$seqid, prodigal_gff$ID)
setnames(scaffold_to_geneid, old=c("V1", "V2"), new=c("scaffold_id", "gene_id"))
scaffold_to_geneid$gene_no <- tstrsplit(scaffold_to_geneid$gene_id, "_", keep=c(2))
scaffold_to_geneid$prodigal_nt_id <- data.table(paste(scaffold_to_geneid$scaffold_id,"_",scaffold_to_geneid$gene_no, sep=""))
blast_gene_ids<- merge(min_evalues, scaffold_to_geneid, by="prodigal_nt_id", all=TRUE)
blast_gff <- merge(blast_gene_ids, gene_coords, by.x="gene_id", by.y="ID", all=TRUE)
blast_gff_table <- blast_gff[,c(16,17,2,18,19,20,4,5,12,13,14,15,21,22)]
blast_gff_table$annotation <- tstrsplit(blast_gff_table$annotation, "<>", keep=c(1))
sum(!is.na(blast_gff_table$annotation))
length(unique(blast_gff_table$scaffold_id))

final_blast_gff_table <- merge(blast_gff_table, no_hits, by="prodigal_nt_id")

fwrite(blast_gff_table, snakemake@output[['blastp_gff']])
fwrite(min_evalues, snakemake@output[['prodigal_blastp_res_table']])

#write log
sessionInfo()