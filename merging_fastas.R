library(seqinr)
library(msa)  ##this library is the alignment software
library(ape)  # this aids in creating a tree

## make sure to set the working directory to the folder with all the sequences
setwd("~/programs/R/merging_fasta/outbreak1")
############# works by reading in single files but you have to create each vari###############
fa1 = read.fasta(file = "gisaid_hcov-19_2021_08_10_20.fasta", as.string = TRUE)
fa2 = read.fasta(file = "gisaid_hcov-19_2021_08_10_20_a.fasta", as.string = TRUE)


##################### TRYING TO AUTOMATE IT WITH LIST OF FILES ############
files <- list.files(pattern = "\\.fasta$")

DF <-  read.fasta(files[1], as.string = TRUE)

for (f in files[-1]){
  df <- read.fasta((f), as.string = TRUE)      # read the file
  DF <- rbind(DF, df)    # append the current file
  write.fasta(DF, names = getName(DF), file.out = 'merged.fasta')
}

## Reading the document back into R and Checking the names of the sequences
merged_seq = read.fasta(file = "merged.fasta")
getName(merged_seq)


## this step preps the file with all the sequences for the alignment

my_seq = readDNAStringSet("merged.fasta")

#this step creates the alignment.  This step with take a bit of time.
# the step is using the default setting within MSA (clustW)
COValign = msa(my_seq)  
#msaPrettyPrint(COValign, output = 'asis')
# this aids in prepping the document for creating the tree (creates the matrix) 
COValign2 = msaConvert(COValign, type = "seqinr::alignment")

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

alignment2Fasta(COValign, 'outalign.fasta')

d = dist.alignment(COValign2, "identity")
as.matrix(d)

### this is a neighbor joining tree ####

COVtree_nj = nj(d)
plot(COVtree_nj, main="Outbreak #2")
write.tree(COVtree_nj, file = "COVTREE.nwk")


### Minimum Spanning Tree ###
?mst
COVtree_mst = mst(d)
plot(COVtree_mst, main = "Outbreak #2")
write.tree(COVtree_mst, file = "COVTREE_MST.nwk")



###############  playing around ######
add.scale.bar(cex = 0.7, font = 2, col = "red")
layout(1)
plot(COVtree_nj, main="Outbreak #2", type = "unrooted")
plot(COVtree, main="Outbreak #2", type = "fan")
plot(COVtree, main="Outbreak #2", type = "cladogram")
COVtree$tip.label
COVtree$edge
COVtree$Nnode
