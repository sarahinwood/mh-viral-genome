library(data.table)
library(dplyr)
library(rhmmer)

hmmscan_res <- data.table(read_domtblout("output/hmmer/prodigalPFAM_domtblout.out"))
##filter by individual e-value - 5e-04 (same as contig blast)
hmmscan_sig <- filter(hmmscan_res, domain_ievalue<5e-04)
##how many genes have predicted domains - 23
length(unique(hmmscan_sig$query_name))
##how many unique domain predictions - 31
length(unique(hmmscan_sig$domain_name))

##filter important columns out
hmmscan_sig_table <- hmmscan_sig[,c(4,2,23,20,21,13,12)]
fwrite(hmmscan_sig_table, "output/hmmer/sig_domains.csv")

blast_annots <- fread("output/prodigal_blast/blast_annotation_table2.csv")
blast_hmmer <- merge(blast_annots, hmmscan_sig_table, by.x="Prodigal_gene_ID", by.y="query_name", all=TRUE)
fwrite(blast_hmmer, "output/blast_hmmer_res/blast_hmmer_res.csv")
