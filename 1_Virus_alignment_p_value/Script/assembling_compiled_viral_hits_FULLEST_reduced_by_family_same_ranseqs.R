library("BiocManager")
library("BiocGenerics")
library("Biostrings")
library("GenomicRanges")
library("tidyverse")
library("splitstackshape")
library("pryr")
library("stringi")
library("dplyr")
library("devtools")
library("ggpubr")
library("stringr")

setwd("/Users/Martin/Desktop/UROP_WORK/My_code")

cNRE <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/ECRRN_full_code/ECRRN/data/raw/dna_seqs/ECRRN.fasta")

output_fasta_vector_viruses <- vector(mode = "character")

viruses_df <- read.csv("/Users/Martin/Desktop/UROP_WORK/My_code/virus_top_2000_by_alignment_score_same_ranseq_set_10000_ranseqs.csv")

output_fasta_vector_viruses[1] <- as.character(">cNRE")
output_fasta_vector_viruses[2] <- as.character(cNRE)

for (i in seq(from = 3, to = (2*length(viruses_df$virus.name) + 1), by = 2)) {
  output_fasta_vector_viruses[i] <- paste(">", as.character(viruses_df$virus.name[(i - 1)/2]), sep = "")
  output_fasta_vector_viruses[i + 1] <- as.character(viruses_df$virus.alignment.hit.sequence[(i - 1)/2])
  print(i)
  print(i + 1)
}

writeLines(output_fasta_vector_viruses, "top_sequences_by_family_same_ranseqs.txt", sep = "\n", useBytes = FALSE)



