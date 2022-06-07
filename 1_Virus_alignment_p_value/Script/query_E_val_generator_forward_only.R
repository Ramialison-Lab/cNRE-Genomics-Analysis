query_E_val_generator_forward_only <- function (query, genome_database) {
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
  library("ggplot2")
  
  #query_reversecomp <- reverseComplement(query)
  
  query_character <- as.character(query)
  #query_reversecomp_character <- as.character(query_reversecomp)
  
  query_dnastringset <- DNAStringSet(query_character)
  #query_reversecomp_dnastringset <- DNAStringSet(query_reversecomp_character)
  
  genome_character <- as.character(genome_database)
  genome_dnastringset <- DNAStringSet(genome_character)
  
  query_score <- score(pairwiseAlignment(genome_dnastringset, query_dnastringset, type = "local-global"))
  #query_reversecomp_score <- score(pairwiseAlignment(genome_dnastringset, query_reversecomp_dnastringset, type = "local-global"))
  
  random_sequence_pairwiseAlignment <- vector(mode = "list", length = 1000)
  random_sequence_score <- vector(length = 1000)
  
  for (m in 1:1000) {
    temp_ranseq <- stri_rand_strings(1, 32, pattern = "[A T C G]")
    
    temp_ranseq_c <- as.character(temp_ranseq)
    temp_ranseq_S4 <- DNAStringSet(temp_ranseq_c)
    ranseq_pwisealignment <- pairwiseAlignment(genome_dnastringset, temp_ranseq_S4, type = "local-global")
    ranseq_alscore <- score(ranseq_pwisealignment)
    
    random_sequence_score[m] <- ranseq_alscore
    print(m)
  }
  
  random_sequence_score_desc <- sort(random_sequence_score, decreasing = TRUE)
  
  j = 1
  
  if (query_score >= random_sequence_score_desc[j]) {
    p_forward = j/length(random_sequence_score_desc)
    print(j)
  } else {
    while(query_score < random_sequence_score_desc[j]) { #This is the line with the error in it!
      j = j + 1
      print(j)
    }
    p_forward = j/length(random_sequence_score_desc)
  }
  
  #if (query_reversecomp_score >= random_sequence_score_desc[k]) {
    #p_reverse = k/length(random_sequence_score_desc)
    #print(k)
  #} else {
    #while(query_reversecomp_score < random_sequence_score_desc[k]) {
      #k = k + 1
      #print(k)
    #}
   # p_reverse = k/length(random_sequence_score_desc)
  #}
  
  #output_list <- c(p_forward, p_reverse)
  
  return(p_forward)
  
}










 