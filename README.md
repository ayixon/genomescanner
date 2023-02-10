${\color{White}Fast  \space genome \space classifier  \space is \space a  \space tool \space that \space can  \space quickly \space  and  \space accurately \space classify \space prokaryotic  \space genomes}$


  Experimental tool to taxonomically classify microbial genomes in a few minutes

COPYRIGHT="Copyright (C) 2022 Ayixon Sánchez Reyes"                                                               
This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU         
General Public License as published by the Free Software Foundation.This program is distributed WITHOUT ANY WARRANTY.                                              

```diff
+DEPENDENCIES: Mash; JolyTree, apcalc, ncbi-entrez-direct, orthoani, Blast, Biopython                               
```                                                                                                               

**Before you begin, install the following**:                                                                 

    sudo apt install mash

    sudo apt install apcalc

    sudo apt install ncbi-entrez-direct

    sudo apt install ncbi-blast+ 

    pip install orthoani

    conda install -c bioconda jolytree
         
  **Usage**:  
  
     ./GenoScanner.sh -i query_genome	
     
     
```diff
- The working directory must contain the mash database (.msh) and the query genome in fasta format

```
> **Download preconfigured** database here [Mash DB format .msh](https://figshare.com/ndownloader/files/37939296)

  > This database contains ~18,000 genomic records with standing in nomenclature

Rational: Compare a query_genome vs a curated MASH database;  select the nearest phylogenetic neighbors; 
Estimate the ANI of the query vs the references; store the genomes in a folder and pass them to JolyTree for phylogenetic estimation.                                   

Fast genome classifier deals with the "Phylophenetic Species Concept" by testing two of its hypotheses:

    The Genomic Coherence measured through the genomic distance of Mash and ANI
     
    The phylogenetic hypothesis of monophyly
        
```diff
@@ Ayixon Sánchez-Reyes   ayixon@gmail.com @@

Computational Microbiology   

Microbiological Observatory 

Institute of Biotechnology, UNAM, Cuernavaca, MEXICO 

VERSION=1.0. Written by Ayixon Sánchez Reyes    
###########################################################################################################               
