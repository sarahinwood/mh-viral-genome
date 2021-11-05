library(data.table)
library(dplyr)
library(rtracklayer)

prodigal_gff <- readGFF("output/prodigal/spades/gene_predictions.gff")
prodigal_blastp <- fread("output/prodigal/spades/prodigal_blastp.outfmt6")

#######################
## prodgal gff file  ##
#######################

prodigal_gff$conf <-as.numeric(prodigal_gff$conf)
setorder(prodigal_gff, -conf)

##complete genes
complete <- sum(prodigal_gff$partial=="00")
cat(complete, "complete genes")
##incomplete at right edge
incomplete_r <- sum(prodigal_gff$partial=="10")
cat(incomplete_r, "incomplete at right-edge genes")
##incomplete at left edge
incomplete_l <- sum(prodigal_gff$partial=="01")
cat(incomplete_l, "incomplete at left-edge genes")
##incomplete both edges
incomplete_rl <- sum(prodigal_gff$partial=="11")
cat(incomplete_rl, "incomplete at both edge genes")

##conf: A confidence score for this gene, representing the probability that this gene is real, i.e. 78.3% means Prodigal believes that gene is real 78.3% of the time and a false positive 21.7% of the time.
##score: The total score for this gene.
##cscore: The hexamer coding portion of the score, i.e. how much this gene looks like a true protein.
##sscore: A score for the translation initiation site for this gene; it is the sum of the following three fields.
##rscore: A score for the RBS motif of this gene.
##uscore: A score for the sequence surrounding the start codon.
##tscore: A score for the start codon type (ATG vs. GTG vs. TTG vs. Nonstandard).
##mscore: A score for the remaining signals (stop codon type and leading/lagging strand information).

################
## blastp res ##
################

setnames(prodigal_blastp, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14"),
         new=c("prodigal_nt_id", "nr_db_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings",
               "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "taxid", "annotation"))
##order so that in event of eval min. tie, which.min takes hit with highest bitscore
setorder(prodigal_blastp, prodigal_nt_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide, sorted by bit-value in case of e-value tie
min_evalues <- prodigal_blastp[,.SD[which.min(evalue)], by=prodigal_nt_id]

####################################
## merge blastp res with prodigal ##
####################################

gene_coords <- subset(prodigal_gff, select=c(start, end, ID, conf, partial))
##scaffold id to gene id
scaffold_to_geneid <- data.table(prodigal_gff$seqid, prodigal_gff$ID)
setnames(scaffold_to_geneid, old=c("V1", "V2"), new=c("scaffold_id", "gene_id"))
scaffold_to_geneid$gene_no <- tstrsplit(scaffold_to_geneid$gene_id, "_", keep=c(2))
scaffold_to_geneid$prodigal_nt_id <- data.table(paste(scaffold_to_geneid$scaffold_id,"_",scaffold_to_geneid$gene_no, sep=""))
blast_gene_ids<- merge(min_evalues, scaffold_to_geneid, by="prodigal_nt_id", all=TRUE)
blast_gff <- merge(blast_gene_ids, gene_coords, by.x="gene_id", by.y="ID", all=TRUE)
blast_gff_table <- blast_gff[,c(16,17,2,18,19,20,21,4,12,13,14,15)]
blast_gff_table$annotation <- tstrsplit(blast_gff_table$annotation, "<>", keep=c(1))
blast_gff_table$scaffold_no <- tstrsplit(blast_gff_table$scaffold_id, "_", keep=c(2))

fwrite(blast_gff_table, snakemake@output[['blastp_gff']])
fwrite(min_evalues, snakemake@output[['prodigal_blastp_res_table']])