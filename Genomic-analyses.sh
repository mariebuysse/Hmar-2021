# NCBI COG pipeline
# We use the proteomes filtered from Pseudofinder, using the <GENOME>_intact.
# Each proteome should be renamed as <GENOME>.faa and the gene inside as <GENOME>_<NUMBER> e.g. FLEESP_1 gene inside FLEESP.faa 
# We should be in a folder with all the .faa to analyse together with COG pipeline
cat *.faa >  master_genomes
grep "^>" master_genomes| perl -pe 's/^>([A-Za-z]+\d*)(_\d+)/\1\2,\1/' > GenQuery.p2o.csv
cat GenQuery.p2o.csv cog2003-2014.csv > tmp.p2o.csv
mkdir BLASTss BLASTno BLASTff BLASTcogn
makeblastdb -dbtype prot -in master_genomes
psiblast -query master_genomes -db master_genomes -show_gis -outfmt 7 -num_alignments 10 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTss/QuerySelf.tab
psiblast -query master_genomes -db <COG_db> -show_gis -outfmt 7 -num_alignments 1000 -dbsize 100000000 -comp_based_stats F -seg no -out BLASTno/QueryCOGs.tab
psiblast -query master_genomes -db <COG_db> -show_gis -outfmt 7 -num_alignments 1000 -dbsize 100000000 -comp_based_stats T -seg yes -out BLASTff/QueryCOGs.tab
COG/COGmakehash/COGmakehash -i=tmp.p2o.csv -o=./BLASTcogn -s="," -n=1
COG/COGreadblast/COGreadblast -d=./BLASTcogn -u=./BLASTno -f=./BLASTff -s=./BLASTss -e=0.1 -q=2 -t=2
COG/COGcognitor/COGcognitor -i=./BLASTcogn -t=cog2003-2014.csv -q=GenQuery.p2o.csv -o=GenQuery.COG.csv


# Launch orthofinder, using blast
orthofinder -f <genomes_location/> -S blast
# Go to the results folder of Orthofinder
cd Orthofinder/<Results>/
# We get a subset of the table with the orthogroups and the genome name (It requires " Prokka-styled" gene names, i.e. "GENOME_$")
tail -n +2 Phylogenetic_Hierarchical_Orthogroups/N0.tsv | awk '{ $1=""; $3=""; print $0 }' | sed 's/,//g' | perl -pe 's/([A-Z]+)_[0-9]+/\1/gi' > orthogroups.tsv
# We add the genes unassigned to orthogroups
tail -n +2 ../Orthogroups/Orthogroups_UnassignedGenes.tsv | sed 's/,//g' | perl -pe 's/([A-Z]+)_[0-9]+/\1/gi' >> orthogroups.tsv
# We use the table to extract the list of orthogroups present for each genome
python3 splitta_cogtab.py orthogroups.tsv <GENOME> | sort | uniq > <GENOME>.og

# script "splitta_cogtab.sh"
# name to search from table
# use as: python3 splitta_cogtab.py FLEISR orthogroups.tsv
import sys
import re
genome = sys.argv[2]
orthofile = sys.argv[1]
f = open(orthofile)
for line in f:
	linea=line.strip().split("\t")
	for l in linea[1:]:
		if l == genome:
			print (linea[0])


# Pseudogene prediction, using Prokka, Pseudofinder, and hoc script to get a subset of genes
# These steps were repeated with each genome
# We launched Prokka on each genome to analyse to generate a compliant genbank file, required by pseudofinder
prokka --rfam --compliant <FASTA>
# We launched pseudofinder 1.0 on each genome
python3 /pseudofinder-1.0/pseudofinder.py annotate --genome <genbank file> --outprefix <prefix> --diamond -db <path to diamond db>
# $1 sigla tipo sigla_intact.faa
# $2 ffn
grep '^>' $1_intact.faa | sed 's/>//g' | awk '{print $1}' > lista_buoni
python3 ~/script/acchiappa_fasta.py $1_proteome.faa lista_buoni -r	
grep '^>' $1_proteome.faa.subset.fasta | sed 's/>//g' | awk '{print $1}' > lista_scarti
python3 ~/script/acchiappa_fasta.py $2 lista_scarti -r
# Extract the list of intact proteins from the outputs of pseudofinder
grep '^>' <prefix>_intact.faa | sed 's/>//g' | awk '{print $1}' > genes2keep.list
# The ad hoc script will take all genes not present in the list
python3 get_fasta.py <prefix>_proteome.faa genes2keep.list -r
# Using the proteome as target we got "<prefix>_proteome.faa.subset.fasta" , a fasta with the proteins predicted as pseudogenes
# We extract the gene names from the fasta
grep '^>' <prefix>_proteome.faa.subset.fasta | sed 's/>//g' | awk '{print $1}' > genes2rumenta.list
# We use this list to get a subset from <genome.ffn> made by Prokka
python3 get_fasta.py <genome>.ffn genes2rumenta.list -r
# <genome>.ffn.subset.fasta is the fasta containing the DNA sequences of all genes, including RNAs, predicted as intact

# script "get-fasta.sh"
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:32:57 2019
"""
from Bio import SeqIO
import sys
with open(sys.argv[2]) as c:
	lista_nomi=c.read().splitlines()
out = (str(sys.argv[1])+".subset.fasta")
if len(sys.argv)==4:
		if sys.argv[3]=="-r":
			print("Excluding the list")
			records = (r for r in SeqIO.parse(sys.argv[1], "fasta") if r.id not in lista_nomi)
		else:
			print("Wrong sysargv")
else:
	records = (r for r in SeqIO.parse(sys.argv[1], "fasta") if r.id in lista_nomi)
SeqIO.write(records, out, "fasta")


# pseudogenization level representation (python)
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import pandas as pd
import sys
data=sys.argv[1]
tab = pd.read_csv(data, index_col=0,sep="\t").T
tab_minmax=(tab/tab.max())
print(tab_minmax)
colour = sns.color_palette("rocket_r", as_cmap=True)
sns.heatmap(tab_minmax, robust=True, xticklabels=1,cmap=colour,linewidths=.5)
plt.show()
