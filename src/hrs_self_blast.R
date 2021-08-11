library(data.table)
library(dplyr)

self_blast <- fread("output/hrs_self_blast/blastn.outfmt6")
contig_lengths <- fread("data/viral_contig_lengths.csv")

setnames(self_blast, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

##will contain each hit twice though
nonself_hits <- subset(self_blast, !(self_blast$query==self_blast$subject))
##merge with table of contig lengths to see whether near edges? or whether in middle?
hits_lengths <- merge(nonself_hits, contig_lengths, by.x="query", by.y="contig_id", all.x=TRUE)
hits_lengths_table <- hits_lengths[,c(1,2,3,4,14,7,8,9,10,11)]

##don't pick up any repeat content on edges of contigs as they did in LbFV
##all regions ID'd by blast are in regions with predicted ORFs that have high homology to Bro genes
##LbFV only 1 Bro gene, most viruses have multiple though (LbFV ref?)

##how many predictions contain strings of Xs
##how many contigs contain strings of Ns
