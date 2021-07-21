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
python3 split_cogtab.py orthogroups.tsv <GENOME> | sort | uniq > <GENOME>.og

# Pseudogene prediction, using Prokka, Pseudofinder, and hoc script to get a subset of genes
# These steps were repeated with each genome
# We launched Prokka on each genome to analyse to generate a compliant genbank file, required by pseudofinder
prokka --rfam --compliant <FASTA>
# We launched pseudofinder 1.0 on each genome
python3 /pseudofinder-1.0/pseudofinder.py annotate --genome <genbank file> --outprefix <prefix> --diamond -db <path to diamond db>
# $1 prefix
# $2 ffn
grep '^>' $1_intact.faa | sed 's/>//g' | awk '{print $1}' > lista_buoni
python3 get_fasta.py  $1_proteome.faa lista_buoni -r	
grep '^>' $1_proteome.faa.subset.fasta | sed 's/>//g' | awk '{print $1}' > lista_scarti
python3 get_fasta.py  $2 lista_scarti -r
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

##FB vitamins, heme, FPI pathways and mismatch repair system analyses
makeblastdb -in finalgenome_scaffolds.fasta -dbtype nucl -out genome_db
blastn -query query-gene.fasta -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -perc_identity 40 -out gene_vs_genome.out -db genome_db 
makeblastdb -in annotation.faa -dbtype prot -out genome_db
blastp -query query-gene.faa -out gene_vs_genome.out -num_threads 6 -db genome_db -outfmt "6 qseqid sseqid sseq qlen pident nident mismatch evalue sstart send gapopen" -evalue 1e-1
