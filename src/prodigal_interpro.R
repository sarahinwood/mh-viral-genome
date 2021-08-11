#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

interpro <- fread(snakemake@input[['interpro_res']], fill=TRUE)

setnames(interpro, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("prodigal_nt_id", "MD5_digest", "seq_length", "analysis", "sig_accession", "sig_description",
               "start", "stop", "evalue", "status", "date", "interpro_accession", "interpro_description", "GO"))
interpro_table <- interpro[,c(1,6,9,13,14)]
fwrite(interpro_table, snakemake@output[['interpro_table']])