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

setwd("/Users/Martin/Desktop/UROP_WORK/My_code")

output_fasta_vector <- vector(mode = "character")

Taeniopygia_guttata <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/taeniopygia_guttata_genome.fna")
temp_taeGut <-reverseComplement(Taeniopygia_guttata$`CM000515.1 Taeniopygia guttata chromosome 1, whole genome shotgun sequence`[90126026:90126057])
length(temp_taeGut)
temp_taeGut_character <- as.character(temp_taeGut)
output_fasta_vector[1] <- "> Taeniopygia_guttata"
output_fasta_vector[2] <- temp_taeGut_character
rm(Taeniopygia_guttata)

Numida_meleagris <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/numida_meleagris_genome.fna")
temp_numMel <- Numida_meleagris$`CM007831.1 Numida meleagris isolate 19003 breed g44 Domestic line chromosome 18, whole genome shotgun sequence`[23819:23850]
length(temp_numMel)
temp_numMel_character <- as.character(temp_numMel)
output_fasta_vector[3] <- "> Numida_meleagris"
output_fasta_vector[4] <- temp_numMel_character
rm(Numida_meleagris)

Tympanuchus_cupido <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/tympanuchus_cupido_pinnatus_genome.fna")
temp_tymCup <-reverseComplement(Tympanuchus_cupido$`MOXI01000045.1 Tympanuchus cupido pinnatus isolate GPC 3440 Scf7LSr2_45, whole genome shotgun sequence`[5866619:5866650])
length(temp_tymCup)
temp_tymCup_character <- as.character(temp_tymCup)
output_fasta_vector[5] <- "> Tympanuchus_cupido"
output_fasta_vector[6] <- temp_tymCup_character
rm(Tympanuchus_cupido)

Lyrurus_tetrix_tetrix <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/lyrurus_tetrix_tetrix_genome.fna")
temp_lyrTetTet <- Lyrurus_tetrix_tetrix$`JDSL01324135.1 Lyrurus tetrix tetrix chr19contig102, whole genome shotgun sequence`[151:182]
length(temp_lyrTetTet)
temp_lyrTetTet_character <- as.character(temp_lyrTetTet)
output_fasta_vector[7] <- "> Lyrurus_tetrix_tetrix"
output_fasta_vector[8] <- temp_lyrTetTet_character
rm(Lyrurus_tetrix_tetrix)

Meleagris_gallopavo <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/meleagris_gallopavo_genome.fna")
temp_melGal <- reverseComplement(Meleagris_gallopavo$`KN389535.1 Meleagris gallopavo isolate NT-WF06-2002-E0010 breed Aviagen turkey brand Nicholas breeding stock unplaced genomic scaffold ChrUn_random_7180001957605, whole genome shotgun sequence`[15013:15044])
length(temp_melGal)
temp_melGal_character <- as.character(temp_melGal)
output_fasta_vector[9] <- "> Meleagris_gallopavo"
output_fasta_vector[10] <- temp_melGal_character
rm(Meleagris_gallopavo)

Gallus_gallus <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/gallus_gallus.fna")
temp_galGal <- reverseComplement(Gallus_gallus$`NC_006106.5 Gallus gallus breed Red Jungle Fowl isolate RJF #256 chromosome 19, GRCg6a`[72671:72702])
length(temp_galGal)
temp_galGal_character <- as.character(temp_galGal)
output_fasta_vector[11] <- "> Gallus_gallus"
output_fasta_vector[12] <- temp_galGal_character
rm(Gallus_gallus)

Pavo_cristatus <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Pavo_cristatus.fna")
temp_pavCri <- reverseComplement(Pavo_cristatus$`QZWQ01000145.1 Pavo cristatus isolate SKPea2016_SI scaffold202_len342489_cov0, whole genome shotgun sequence`[337180:337211])
length(temp_pavCri)
temp_pavCri_character <- as.character(temp_pavCri)
output_fasta_vector[13] <- "> Pavo_cristatus"
output_fasta_vector[14] <- temp_pavCri_character
rm(Pavo_cristatus)

Chrysolophus_pictus <-readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Chrysolophus_pictus.fna")
temp_chrPic <- reverseComplement(Chrysolophus_pictus$`JRFT01060738.1 Chrysolophus pictus breed wild golden pheasant contig60738, whole genome shotgun sequence`[904:935])
length(temp_chrPic)
temp_chrPic_character <- as.character(temp_chrPic)
output_fasta_vector[15] <- "> Chrysolophus_pictus"
output_fasta_vector[16] <- temp_chrPic_character
rm(Chrysolophus_pictus)

Phasianus_colchicus <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Phasianus_colchicus.fna")
temp_phaCol <- reverseComplement(Phasianus_colchicus$`QCWP01000047.1 Phasianus colchicus isolate SZU-A5-319 scaffold55_cov57, whole genome shotgun sequence`[9855508:9855539])
length(temp_phaCol)
temp_phaCol_character <- as.character(temp_phaCol)
output_fasta_vector[17] <- "> Phasianus_colchicus"
output_fasta_vector[18] <- temp_phaCol_character
rm(Phasianus_colchicus)

Syrmaticus_mikado <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Syrmaticus_mikado.fna")
temp_syrMik <-reverseComplement(Syrmaticus_mikado$`QGNR01000125.1 Syrmaticus mikado isolate TWDZ01 scaffold125, whole genome shotgun sequence`[24073:24104])
length(temp_syrMik)
temp_syrMik_character <- as.character(temp_syrMik)
output_fasta_vector[19] <- "> Syrmaticus_mikado"
output_fasta_vector[20] <- temp_syrMik_character
rm(Syrmaticus_mikado)

Coturnix_japonica <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/coturnix_japonica_genome.fna")
temp_cotJap <- reverseComplement(Coturnix_japonica$`KQ966789.1 Coturnix japonica isolate 7356 unplaced genomic scaffold chrUnrandom599, whole genome shotgun sequence`[13552:13581])
length(temp_cotJap)
temp_cotJap_character <- as.character(temp_cotJap)
output_fasta_vector[21] <- "> Coturnix_japonica"
output_fasta_vector[22] <- temp_cotJap_character
rm(Coturnix_japonica)

Bambusicola_thoracicus <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Bambusicola_thoracicus.fna")
temp_bamTho <- reverseComplement(Bambusicola_thoracicus$`PPHD01031036.1 Bambusicola thoracicus isolate RTK389 B.A.C_31035, whole genome shotgun sequence`[9246:9277])
length(temp_bamTho)
temp_bamTho_character <- as.character(temp_bamTho)
output_fasta_vector[23] <- "> Bambusicola_thoracicus"
output_fasta_vector[24] <- temp_bamTho_character
rm(Bambusicola_thoracicus)

Callipepla_squamata <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/callipepla_squamata_genome.fna")
temp_calSqu <- Callipepla_squamata$`MCFN01001451.1 Callipepla squamata strain Texas jcf7180005229139, whole genome shotgun sequence`[12433:12464]
length(temp_calSqu)
temp_calSqu_character <- as.character(temp_calSqu)
output_fasta_vector[25] <- "> Callipepla_squamata"
output_fasta_vector[26] <- temp_calSqu_character
rm(Callipepla_squamata)

Colinus_virginianus <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/colinus_virginianus_genome.fna")
temp_colVir <- reverseComplement(Colinus_virginianus$`AWGT02003190.1 Colinus virginianus strain Texas jcf7180006208378, whole genome shotgun sequence`[11965:11996])
length(temp_colVir)
temp_colVir_character <- as.character(temp_colVir)
output_fasta_vector[27] <- "> Colinus_virginianus"
output_fasta_vector[28] <- temp_colVir_character
rm(Colinus_virginianus)

Anas_platyrhynchos <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Anas_platyrhyncho.fna")
temp_anaPla <- reverseComplement(Anas_platyrhynchos$`CM011859.1 Anas platyrhynchos isolate PK-2015 chromosome 26, whole genome shotgun sequence`[1283767:1283798])
length(temp_anaPla)
temp_anaPla_character <- as.character(temp_anaPla)
output_fasta_vector[29] <- "> Anas_platyrhynchos"
output_fasta_vector[30] <- temp_anaPla_character
rm(Anas_platyrhynchos)

Anas_zonorhyncha <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Anas_zonorhyncha.fna")
temp_anaZon <- Anas_zonorhyncha$`NOIK01000635.1 Anas zonorhyncha breed spot-billed scaffold637, whole genome shotgun sequence`[268033:268064]
length(temp_anaZon)
temp_anaZon_character <- as.character(temp_anaZon)
output_fasta_vector[31] <- "> Anas_zonorhyncha"
output_fasta_vector[32] <- temp_anaZon_character
rm(Anas_zonorhyncha)

Cairina_moschata <- readDNAStringSet("/Users/Martin/Desktop/UROP_WORK/BLAST_analysis_data/genomes/further_genomes/Cairina_moschata_domestica.fna")
temp_caiMos <- reverseComplement(Cairina_moschata$`QZEJ01000539.1 Cairina moschata domestica isolate CanardBarbarie3_02_2015_G Scaffold539, whole genome shotgun sequence`[104937:104968])
length(temp_caiMos)
temp_caiMos_character <- as.character(temp_caiMos)
output_fasta_vector[33] <- "> Cairina_moschata"
output_fasta_vector[34] <- temp_caiMos_character
rm(Cairina_moschata)

writeLines(output_fasta_vector, "avian_cNRE_alignments.txt", sep = "\n", useBytes = FALSE)
#This file can be converted by hand into a fasta file by simply changing the .txt to a .fasta
#The cNRE was then manually added to the file (under Coturnix Japonica)
