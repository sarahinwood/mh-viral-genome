library(rtracklayer)
library(data.table)

#Mh viral contigs bbstats length
mh_viral_contig_length <- 92650

prodigal_gff <- data.frame(readGFF("data/Mh_prodigal/gene_predictions.gff"))

##nt length ranges 108-3867 nts, mean nt length =  866.6bp
prodigal_gff$nt_length <- (prodigal_gff$end)-(prodigal_gff$start)+1
mean(prodigal_gff$nt_length)

## AA length range 36-1289, mean AA length = 288.9 aas
prodigal_gff$aa_length <- prodigal_gff$nt_length/3
mean(prodigal_gff$aa_length)

#91.7% coding density in 92.6kb total (LbFV 80.0% in 111kb)
sum(prodigal_gff$nt_length)/mh_viral_contig_length

#98 genes, 22 incomplete at one or both edges though)
sum(!(prodigal_gff$partial=="00"))
#5 genes are complete but don't start with ATG - viruses only ATG

#33.8% GC content, 66.2% AT (LbFV 78.7% AT)
#GC content of each contig?
#depth of each contig?

#50 genes on - strand
neg_strand <- subset(prodigal_gff, strand=="-")
#total length = 40.4kb
sum(neg_strand$nt_length)

#48 genes on + strand
pos_strand <- subset(prodigal_gff, strand=="+")
#total length = 44.5kb
sum(pos_strand$nt_length)

##no of genes with bro homology - LbFV has one - what does this mean?

##repeat content? on edges of contigs?

#how many of these genes have blast hits
#how many have LbFV blast hits

##spread of genes on either strand
