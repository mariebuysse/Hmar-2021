orthofinder -f ./symbiont_phylo -t 4 -S blast ##symbiont_phylo: directory containing all genomes of FLE or Midichloria or Rickettsia used for the phylogenomic analyses 
# ADD command line MUSCLE and Gblocks
modeltest-ng -i symbiont-phylo_concatenated.faa -p 12 -T raxml -d aa
raxmlHPC-PTHREADS -s Selected_AlnGblock.fasta -m model -n Selected_AlnGblock.fasta.boot -x 1234 -p 123 -f a -T 12 -# 100 
