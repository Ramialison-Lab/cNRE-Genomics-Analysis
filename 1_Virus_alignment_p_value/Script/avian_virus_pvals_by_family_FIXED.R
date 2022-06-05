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
source("query_E_val_generator_forward_only.R")
cNRE <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/ECRRN_full_code/ECRRN/data/raw/dna_seqs/ECRRN.fasta")
avian_viruses_full_set_ncbi <- readDNAStringSet("C:/Users/Martin/Desktop/UROP_WORK/My_code/avian_viruses_complete_19_2_20.fasta")

input_df <- read.csv("/Users/Martin/Desktop/UROP_WORK/My_code/compiled_top_2000_avian_viruses_by_family_reduced.csv")
eval_list <- vector(mode = "list", length = nrow(input_df))

ranseq_vector <- vector(mode = "list", length=10000)

for (i in 1:length(ranseq_vector)) {
  temp_ranseq <- stri_rand_strings(1, 32, pattern = "[A T C G]")
  temp_ranseq_c <- as.character(temp_ranseq)
  temp_ranseq_S4 <- DNAStringSet(temp_ranseq_c)
  ranseq_vector[i] = temp_ranseq_S4
  print(i)
}

for (i in 1:nrow(input_df)) {
  n = 1
  if (isTRUE(input_df$strand.type[i] == "reverse complement")) {
    temp_index <- input_df$index.in.avian_viruses_complete_19_2_20.fasta[i]
    temp_sequence <- DNAStringSet(avian_viruses_full_set_ncbi[temp_index], use.names = FALSE)
    #NEW SCRIPT
    temp_alscore_vector = vector(length=length(ranseq_vector))
    for (j in 1:length(ranseq_vector)){
      temp_alscore_vector[j] = score(pairwiseAlignment(temp_sequence, ranseq_vector[[j]], type = "local-global"))
      print(n)
      n = n + 1
    }
    temp_alscore_vector_desc <- sort(temp_alscore_vector, decreasing = TRUE)
    
    j = 1
    
    if (input_df$virus.alignment.score[i] >= temp_alscore_vector_desc[j]) {
      p_forward = j/length(temp_alscore_vector_desc)
      print(j)
    } else {
      while(input_df$virus.alignment.score[i] < temp_alscore_vector_desc[j]) { #This is the line with the error in it!
        j = j + 1
        print(j)
      }
      p_forward = j/length(temp_alscore_vector_desc)
    }
    
    eval_list[i] <- p_forward
    print(paste0('Virus number ', i, ' completed'))
    
  } else {
    temp_index <- input_df$index.in.avian_viruses_complete_19_2_20.fasta[i]
    temp_sequence <- DNAStringSet(avian_viruses_full_set_ncbi[temp_index], use.names = FALSE)
    #NEW SCRIPT
    temp_alscore_vector = vector(length=length(ranseq_vector))
    for (j in 1:length(ranseq_vector)){
      temp_alscore_vector[j] = score(pairwiseAlignment(temp_sequence, ranseq_vector[[j]], type = "local-global"))
      print(n)
      n = n + 1
    }
    temp_alscore_vector_desc <- sort(temp_alscore_vector, decreasing = TRUE)
    
    j = 1
    
    if (input_df$virus.alignment.score[i] >= temp_alscore_vector_desc[j]) {
      p_forward = j/length(temp_alscore_vector_desc)
      print(j)
    } else {
      while(input_df$virus.alignment.score[i] < temp_alscore_vector_desc[j]) { #This is the line with the error in it!
        j = j + 1
        print(j)
      }
      p_forward = j/length(temp_alscore_vector_desc)
    }
    
    eval_list[i] <- p_forward
    print(paste0('Virus number ', i, ' completed'))
  }
}






virus_E_value <- unlist(eval_list)
#eval_list_vector <- unlist(eval_list)

output_df <- data.frame(cbind(input_df, virus_E_value))
output_df_character <- mutate_all(output_df, as.character)

write.csv(output_df, "virus_top_2000_by_alignment_score_same_ranseq_set_10000_ranseqs.csv")

#TESTING FOR LOOP FOR 
i = 1
if (isTRUE(input_df$strand.type[i] == "reverse complement")) {
  n = 1
  temp_index <- input_df$index.in.avian_viruses_complete_19_2_20.fasta[i]
  temp_sequence <- DNAStringSet(avian_viruses_full_set_ncbi[temp_index], use.names = FALSE)
  #NEW SCRIPT
  temp_alscore_vector = vector(length=length(ranseq_vector))
  for (j in 1:length(ranseq_vector)){
    temp_alscore_vector[j] = score(pairwiseAlignment(temp_sequence, ranseq_vector[[j]], type = "local-global"))
    print(n)
    n = n + 1
  }
  temp_alscore_vector_desc <- sort(temp_alscore_vector, decreasing = TRUE)
  
  j = 1
  
  if (input_df$virus.alignment.score[i] >= temp_alscore_vector_desc[j]) {
    p_forward = j/length(temp_alscore_vector_desc)
    print(j)
  } else {
    while(input_df$virus.alignment.score[i] < temp_alscore_vector_desc[j]) { #This is the line with the error in it!
      j = j + 1
      print(j)
    }
    p_forward = j/length(temp_alscore_vector_desc)
  }
  
  eval_list[i] <- p_forward
  print(paste0('Virus number ', i, ' completed'))
  
} else {
  temp_index <- input_df$index.in.avian_viruses_complete_19_2_20.fasta[i]
  temp_sequence <- DNAStringSet(avian_viruses_full_set_ncbi[temp_index], use.names = FALSE)
  #NEW SCRIPT
  temp_alscore_vector = vector(length=length(ranseq_vector))
  for (j in 1:length(ranseq_vector)){
    temp_alscore_vector[j] = score(pairwiseAlignment(temp_sequence, ranseq_vector[[j]], type = "local-global"))
    print(n)
    n = n + 1
  }
  temp_alscore_vector_desc <- sort(temp_alscore_vector, decreasing = TRUE)
  
  j = 1
  
  if (input_df$virus.alignment.score[i] >= temp_alscore_vector_desc[j]) {
    p_forward = j/length(temp_alscore_vector_desc)
    print(j)
  } else {
    while(input_df$virus.alignment.score[i] < temp_alscore_vector_desc[j]) { #This is the line with the error in it!
      j = j + 1
      print(j)
    }
    p_forward = j/length(temp_alscore_vector_desc)
  }
  
  eval_list[i] <- p_forward
  print(paste0('Virus number ', i, ' completed'))
  
  
}

