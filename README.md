# cNRE-Genomics-Analysis
 cNRE comparative genomics pipeline
# Comparative genomics analyses of cNRE (complex nuclear receptor cis-regulatory element)
## Description
Computational analyses of cNRE (complex nuclear receptor cis-regulatory element), a 32-bp sequence with hexanucleotide binding repeats in different organisms. 

Developer: **Martin Nikolov**, Ramialison Laboratory, Australian Regenerative Medicine Institute & Murdoch Children's Research Institute, Australia

1. R scripts for calculating custom p-values involved in viral hits alignment (associated with Supplementary Figure 6-7) ("1_Virus_alignment_p_value" directory).

2. Helper R scripts for extracting FASTA sequence from various genomes (associated with Figure 5) ("2_Avian_alignment_extract_FASTA" directory).

Code for the manuscript: Luana Nunes Santos, Ângela Maria da Souza Costa, Martin Nikolov, Allysson Coelho Sampaio, Frank E. Stockdale, Gang F Wangø, Hozana Andrade Castillo, Mariana Bortoletto Grizante, Stefanie Dudczig, Michelle Vasconcelos, Nadia Rosenthal, Patricia Regina Jusuf, Paulo de Oliveira, Tatiana Guimarães de Freitas Matos, William Nikovits Jr., Michael Schubert, Mirana Ramialison*, José Xavier-Neto. Deep origins of the complex Nuclear Receptor Element (cNRE), a cis-regulatory module of viral origin required for preferential expression in the atrial chamber. Under review, preprint available at https://www.biorxiv.org/content/10.1101/2021.11.18.469087v1.full. (<sup>*</sup>: corresponding author)

Languages: R

Operating systems: Windows, Linux, Mac OSX. 

## Getting the Source Code

To get the source code, please click the "fork" button in the upper-right and then add this repo as an upstream source:

````
$ git clone <your_fork_of_the_repo> ppds
$ cd ppds
$ REPO=https://github.com/Ramialison-Lab/cNRE-Genomics-Analysis.git
$ git remote add upstream $REPO
````

To get new content, use 
````
$ git pull upstream master 
````

