orthofinder -f ./symbiont_phylo -t 4 -S blast ##symbiont_phylo: directory containing all genomes of FLE or Midichloria or Rickettsia used for the phylogenomic analyses 

cp -r ./symbiont_phylo/OrthoFinder/Results*/Single_Copy_Orthologue_Sequences ## create a copy of the results directory containing the single copy panorthologs identified by OrthoFinder

for i in Single_Copy_Orthologue_Sequences/*.fa; do muscle -in $i -out $i.aln; done ## align all the proteins in each cluster of single copy panorthologs
for i in Single_Copy_Orthologue_Sequences/*.aln; do Gblocks $i -t=p; done ## remove hypervariable regions in each alignment

# concatenate all the sequences in a single file: symbiont-phylo_concatenated.faa

modeltest-ng -i symbiont-phylo_concatenated.faa -p 12 -T raxml -d aa
raxmlHPC-PTHREADS -s Selected_AlnGblock.fasta -m model -n Selected_AlnGblock.fasta.boot -x 1234 -p 123 -f a -T 12 -# 100 
# FLE : PROTGAMMAICPREV (CPREV+I+G)
# Midichloria : PROTGAMMAILG (LG+I+G)
# Rickettsia : PROTGAMMAIJTT (JTT+G)
