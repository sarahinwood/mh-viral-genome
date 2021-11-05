library(data.table)
library(dplyr)
library(rtracklayer)

blast_res <- fread("output/blast_c4e30_initial/blastn.outfmt6")
original_gff <- readGFF("data/Mh_prodigal/gene_predictions.gff")

setnames(blast_res, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"),
         new=c("query", "subject", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(blast_res, query, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- blast_res[,.SD[which.min(evalue)], by=query]

##how many original genes have blast hit?
length(unique(min_evalues$subject)) ##95

##get list of original gene ids
gene_coords <- subset(original_gff, select=c(start, end, ID, conf, partial, strand))
##scaffold id to gene id
orig_scaffold_to_geneid <- data.table(original_gff$seqid, original_gff$ID)
setnames(orig_scaffold_to_geneid, old=c("V1", "V2"), new=c("scaffold_id", "gene_id"))
orig_scaffold_to_geneid$gene_no <- tstrsplit(orig_scaffold_to_geneid$gene_id, "_", keep=c(2))
orig_scaffold_to_geneid$prodigal_nt_id <- data.table(paste(orig_scaffold_to_geneid$scaffold_id,"_",orig_scaffold_to_geneid$gene_no, sep=""))

##which genes did't have a blast hit?
no_hit <- subset(orig_scaffold_to_geneid, !(prodigal_nt_id %in% min_evalues$subject))

##95/98 original prodigal genes have a hit in new predictions

##three without a blast hit:
      #scaffold_159_11 - fairly short, no Ns, no annot
      #scaffold_252_1 - very short, no Ns, no annot
      #scaffold_1039_2 - reasonable length, no Ns, pithovirus unchar annot

##one gene (scaffold_61_1) had two hits

