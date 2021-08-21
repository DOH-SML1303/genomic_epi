# Reading, combining, and writing FASTA (sequencing) files

# setwd
setwd("~/R/bioinformatics")

# library(seqinr)
library(seqinr)

# read.fasta
fa1 = read.fasta("~/R/bioinformatics/gisaid_hcov-19_2021_08_17_16.fasta", seqtype = c("DNA"), as.string = TRUE)
fa2 = read.fasta("~/R/bioinformatics/gisaid_hcov-19_2021_08_17_16-1.fasta", seqtype = c("DNA"), as.string = TRUE)
fa3 = read.fasta("~/R/bioinformatics/gisaid_hcov-19_2021_08_17_16-2.fasta", seqtype = c("DNA"), as.string = TRUE)

# rbind
fa <- rbind(fa1, fa2, fa3)

# names
names <- as.list(c("seq01", "seq02", "seq03"))

# write.fasta
write.fasta(sequences=fa, names, file.out="sequences.fasta", open = "w")