#!/bin/bash 
	
#FAST GENOME CLASSIFIER IS A TOOL FOR QUICKLY AND ACCURATELY CLASSIFYING PROKARYOTIC GENOMES
####################################################################################################################

####################################################################################################################
#                                                                                                                  #
#COPYRIGHT="Copyright (C) 2022 Ayixon Sánchez Reyes"                                                               #
#                                                                                                                  #
# This program  is free software:  you can  redistribute it  and/or modify it  under the terms  of the GNU         #
# General Public License as published by the Free Software Foundation, either version 3 of the License, or         #
# (at your option) any later version.This program is distributed WITHOUT ANY WARRANTY.                                              #
#                                                                                                                  #
####################################################################################################################
####################################################################################################################
# DEPENDENCIES: Mash; JolyTree, apcalc, ncbi-entrez-direct, orthoani, Blast, Biopython                             #  
#                                                                                                                  #
##########Before you begin, install the following:                                                                 #
#
##########  sudo apt install mash
##########	sudo apt install apcalc
##########	sudo apt install ncbi-entrez-direct
##########  sudo apt install ncbi-blast+ 
########## 	pip install orthoani
########## 	conda install -c bioconda jolytree
#         
####################################################################################################################
#                                                                                                                  #
####################################################################################################################
# Usage: ./GenoScanner.sh -i query_genome	                                                                                   #
# The working directory must contain the mash database (.msh) and the query genome in fasta format                 #
####################################################################################################################

# Rational: Compare a query_genome vs a curated MASH database, select the nearest phylogenetic neighbors;   #
# Estimate the ANI of the query vs the references, store the genomes in a folder                            #
# and pass them to JolyTree for phylogenetic estimation.                                                    #
                                                                                                            #
# FAST GENOME CLASSIFIER deals with the "Phylophenetic Species Concept" by testing two of its hypotheses:   #
     # The Genomic Coherence measured through the genomic distance of Mash and the ANI                      #
     # The phylogenetic hypothesis of monophyly
                                                                                                            #
#  Ayixon Sánchez-Reyes                                                                   ayixon@gmail.com  #
#  Computational Microbiology                                                                               #
#  Microbiological Observatory                                                                              #
#  Institute of Biotechnology, UNAM, Cuernavaca, MEXICO                                                     #
#############################################################################################################
#                                                                                                           #
# ============                                                                                              #
# = VERSIONS =      VERSION=1.0. Written by Ayixon Sánchez Reyes                                            # 
VERSION=1.0                                                                                                 #                       
# ============                                                                                              #
#############################################################################################################
#############################################################################################################

#Mash Distance estimation between the genome query and the database
start=`date +%s.%N`

var2=*.msh 

# Set the variable 'input_file' to an empty string
input_file=''

# Use getopts to parse the input options
while getopts 'i:' opt; do
  case $opt in
    i) input_file=$OPTARG ;;
  esac
done

# Check if the input file is a .fna or .fasta file
if [[ $input_file =~ .*\.(fna|fasta)$ ]]; then
  # Estimate Mash distance with the input file
  mash dist $input_file $var2 -p 10| sort -n -k3 > output.mash.txt; 
echo "" 
else
  # Print an error message if the input file is not a .fna or .fasta file
  echo "Error: Invalid input file. Input file must be a .fna or .fasta file."
fi

echo -e "	\e[0;32mMASH analysis finished. These are the top hits   \e[0m " 

echo "	Query	Reference	D	p_value	shared hashed"

head output.mash.txt


cut -f3 output.mash.txt |head -n 100 > distancias 

uniq distancias > dist.uniq 

echo "" 

echo -e " \e[0;32mPrinting the distance table (Mash D)   \e[0m " 

echo "" 

echo -e "	\e[0;32mThese are the closest genomes to your query and their approximate ANI   \e[0m " 

# Approximate ANI estimation as 1-D
echo "" 

for i in $(fmt dist.uniq)
do  sort -n -k3 output.mash.txt| grep -m1 -F "$i" 
echo "" 
var7=$(calc 1-"$i")
echo -e " \e[0;32mComputing approximate ANI   \e[0m "
echo $var7
echo "---------"  
if ((`bc <<< "$var7>=0.95"`)) ; then
  echo -e '\e[0;33mGenomic Coherence Detected   \e[0m '
  echo "" 
  echo "--------------------------" 
   
fi
  done 
echo "--------" 
echo "" 

#Selection and download of the neighboring genomes with the smallest genomic distance

echo -e " \e[0;32mDownloading the reference genomes from NCBI   \e[0m " 
echo "#####################################" 

for z in $(fmt dist.uniq)
do grep -F -m 1 "$z" output.mash.txt| awk '{print $2}'| sed 's/_[A-Z]/   / gI' | awk '{print $1}' >> Genome_accnumber.txt
  done
   echo ""


cat Genome_accnumber.txt | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget "$url/$fname" ;
        done ;
		
    done
echo "" 
echo -e " \e[0;32mDownload successful   \e[0m "
echo "" 
gzip -d *.gz
#Genomes are copied to a separate directory for phylogenetic analysis
#working directory cleanup operations

echo -e " \e[0;32mStarting Phylogenetic Analysis with JolyTree   \e[0m "
echo "##############################################" 

mkdir Mash_out
cp dist.uniq distancias Mash_out/ 

mkdir JolyTree_in 
cp *.fna JolyTree_in/
cp *.fasta JolyTree_in/
cp *.fa JolyTree_in/

JolyTree.sh -i JolyTree_in/ -b out_tree -t 10 
	
mkdir JolyTree_out
cp out_tree* JolyTree_out/
echo ""
echo -e " \e[0;32m#Phylogenetic analysis finished, results are in the directory JolyTree_out  \e[0m" 
echo ""

echo -e " \e[0;32m#Computing OrthoAni  \e[0m "
echo "--------------------"
echo  Genome_Query: $input_file 

rm *.fna
rm out_tree*

cp Genome_accnumber.txt output.mash.txt Mash_out/ 
rm dist.uniq distancias output.mash.txt Genome_accnumber.txt

cd JolyTree_in/
for ani in *.fna
do orthoani -q $input_file -r $ani -j 4 > $ani.txt
done


awk '{print substr(FILENAME, 3)" \""$0"\" "}' ./*.txt|sort -t\" -nrk2 | sed 's/_[A-Z]/   / gI'| awk '{print $1,"ANI=", $NF}'

cd ..

mkdir OrthoANI_out
cp JolyTree_in/*.txt OrthoANI_out/
echo "" 
rm JolyTree_in/*.txt
rm *.gz
echo ""
	echo -e " \e[0;32mFast_Genome_Classifier pipeline done  \e[0m "
echo ""
echo -e " \e[0;33mWe recommend that you test the phylophenetic hypothesis about *Diagnostic Characters*  \e[0m "

echo -e " \e[0;33mYou may want to explore the following phenotypic determinants:  \e[0m "

echo -e " \e[0;33mMenaquinone biosynthesis including the menaquinone reductase gene (MenJ)  \e[0m "

echo -e " \e[0;33mFatty acid synthesis cycle FASI, FASII \e[0m "

echo -e " \e[0;33m...any other phenotypic determinant of interest  \e[0m "

echo ""
	echo "Thank you for your support, enjoy it"
	echo ""

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc )
echo Time $runtime seconds