library("BiocManager")
library("BiocGenerics")
#BiocManager::install("Biostrings")
library("Biostrings")

rm(list=ls())

cNRE <- readDNAStringSet("ECRRN.fasta")
avian_viruses_full_set_ncbi <- readDNAStringSet("avian_viruses_complete_19_2_20.fasta")
full_promoter <- readDNAStringSet("Full SMyHC III Promoter.fa")

### From Blast source code
##  /** Karlin-Altschul parameter values for substitution scores 1 and -3. */
# static const array_of_8 blastn_values_1_3[] = {
#   { 0, 0, 1.374, 0.711, 1.31, 1.05,  0, 100 },
#   { 2, 2,  1.37,  0.70,  1.2,  1.1,  0,  99 },
#   { 1, 2,  1.35,  0.64,  1.1,  1.2, -1,  98 },
#   { 0, 2,  1.25,  0.42, 0.83,  1.5, -2,  91 },
#   { 2, 1,  1.34,  0.60,  1.1,  1.2, -1,  97 },
#   { 1, 1,  1.21,  0.34, 0.71,  1.7, -2,  88 }  ###Use this row: gap opening = 1, gap extension = 1, lambda=1.21, k=0.34 
# };

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE)
lambda=1.21
k=0.34 
min_p <- 1
min_index <- 0

### Loop through each virus, perform Smith-Waterman alignment and calculate the e-value / probability
for (i in 1:length(avian_viruses_full_set_ncbi)){
  virus <- DNAStringSet(avian_viruses_full_set_ncbi[i], use.names = FALSE)
  S <- pairwiseAlignment(virus, full_promoter, type = "local-global", gapOpening = 1, gapExtension = 1, scoreOnly=T)
  e_value <- k*width(virus)*width(full_promoter)*exp(-lambda*S)
  p_value <- 1 - exp(-e_value)
  if (p_value < min_p){
    min_p <- p_value
    min_index <- i
  }
}






