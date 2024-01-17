```diff
- EXPERIMENTAL tool to taxonomically classify microbial genomes in a few minutes

``` 

COPYRIGHT="Copyright (C) 2022 Ayixon Sánchez Reyes"                                                               
This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU         
General Public License as published by the Free Software Foundation.This program is distributed WITHOUT ANY WARRANTY.                                              

```diff
+DEPENDENCIES: Mash; JolyTree, apcalc, ncbi-entrez-direct, orthoani, Blast, Biopython, bPTP, mptp                                 
```                                                                                                               

**Before you begin, install the following**:                                                                 

    sudo apt update

    sudo apt upgrade -y

    sudo apt install -y build-essential

Install the dependencies

    sudo apt install -y mptp apcalc ncbi-entrez-direct ncbi-blast+

Install Python packages

    pip install orthoani

Install Conda packages

    mamba install -y -c bioconda jolytree

    mamba install -y -c bfurneaux bptp
    
           
  **Usage**:  
  
     ./GenoScanner.sh  -i <input_file> -d <database_file> -m <model>

    "Options:"

    "  -i <input_file>    Input fasta, fna, or fa archive"

    "  -d <database_file> Database file in .msh format"

    "  -m <model>         Select between two models: mptp or bptp"

    "  -h                 Display this help message"
    
```diff
- (bPTP is the default option in our program because it works best for the speciation hypothesis, according to our experience.

```   
   
   
```diff
- The working directory must contain the mash database (.msh) and the query genome in fasta format

```
> **Download preconfigured** database here [Mash DB format .msh](https://figshare.com/ndownloader/files/37939296)

  > This database contains ~18,000 genomic records with standing in nomenclature

Rational: Compare a query_genome vs a curated MASH database;  select the nearest phylogenetic neighbors; 
Estimate the ANI of the query vs the references; store the genomes in a folder and pass them to JolyTree for phylogenetic estimation. Finally, the tree is subjected to speciation hypothesis testing under Poisson Tree Processes Model.                                    

Fast genome classifier deals with the "Phylophenetic Species Concept" by testing following hypotheses:

    The Genomic Coherence measured through the genomic distance of Mash and ANI
     
    The phylogenetic hypothesis of monophyly
    
    The molecular speciation under Poisson Tree Processes Model
        
```diff
@@ Ayixon Sánchez-Reyes   ayixon@gmail.com @@

Computational Microbiology   

Microbiological Observatory 

Institute of Biotechnology, UNAM, Cuernavaca, MEXICO 

VERSION=1.0. Written by Ayixon Sánchez Reyes    
###########################################################################################################               
