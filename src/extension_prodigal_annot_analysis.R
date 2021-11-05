library(tidyverse)
library(data.table)

blast_res_files <- list.files("output/prodigal/price2", pattern="prodigal_blastp_res.csv", recursive=TRUE, full.names = TRUE)
##filenames
names(blast_res_files) <- gsub(".*/.*/.*/(.+)/.*", "\\1", blast_res_files)
##read tables
tables <- lapply(blast_res_files, fread)
names(tables) <- gsub(".*/.*/.*/(.+)/.*", "\\1", blast_res_files)
##full table - separated by filenames
blast_res_table <- rbindlist(tables, idcol="file")

##subset for only viral hits
virus_taxids <- fread("data/species_virus_taxids.txt", header=FALSE)
viral_hits <- subset(blast_res_table, taxid %in% virus_taxids$V1)
viral_counts <- count(viral_hits, file)
setnames(viral_counts, old=c("n"), new=c("viral_counts"))

##subset for non-viral hits
nonviral_hits <- subset(blast_res_table, !(taxid %in% virus_taxids$V1))
nonviral_counts <- count(nonviral_hits, file)
setnames(nonviral_counts, old=c("n"), new=c("nonviral_counts"))

##table of counts for viral and non-viral hits
blast_hit_counts <- merge(viral_counts, nonviral_counts, by="file")

##nonviral hits e4c30
e4.c30.nonviral <- subset(nonviral_hits, file=="extended4.cycle30")


