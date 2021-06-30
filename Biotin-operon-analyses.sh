# phylogeny of the biotin operon
mafft operon-biotin.fas > operon-biotin_aln.fasta
trimal -in operon-biotin_aln.fasta -out operon-biotin_aln_trimal.fasta -fasta -gt 1 -cons 50 
modeltest-ng -i operon-biotin_aln_trimal.fasta -p 12 -T raxml -d aa 
raxmlHPC-PTHREADS -T 8 -f a -s operon-biotin_aln_trimal.fasta -n Phylo_operon.boot -m PROTGAMMAICPREV -x 1234 -p 123 -# 1000 

# bacterial phylogeny
mafft 16S.fas > 16S_aln.fasta
trimal -in 16S_aln.fasta -out 16S_aln_trimal.fasta -fasta -gt 1 -cons 50 
modeltest-ng -i 16S_aln_trimal.fasta -p 12 -T raxml -d aa 
raxmlHPC-PTHREADS -T 8 -f a -s 16S_aln_trimal.fasta -n Phylo_16S.boot -m GTRGAMMAIX -x 1234 -p 123 -# 1000 

## folD and ribH phylogenies
mafft gene.fas > gene_aln.fas
trimal -in gene_aln.fas -out gene_aln_trimal.fasta -fasta -gt 1 -cons 50
modeltest-ng -i gene_aln_trimal.fasta -p 12 -T raxml -d aa 
raxmlHPC-PTHREADS -T 8 -f a -s gene_aln_trimal.fasta -n Phylo_gene.boot -m model -x 1234 -p 123 -# 1000
# model folD: PROTGAMMAICPREV
# model ribH: PROTGAMMALG

# structure of biotin operon (R)
library(ade4)
library(grid)
library(genoPlotR)
library(ggplot2)
symbiontX <- read_dna_seg_from_genbank("symbiontX.gb")

# colors 
# bioA - lightslateblue
# bioD - #FFB740
# bioC - darkslategray1
# bioH - #30C8C4
# bioF - #E57373
# bioB - pink1
symbiontX[2,]$fill[symbiontX[2,]$fill %in% c("blue")]<-"#FFB740"
symbiontX[3,]$fill[symbiontX[3,]$fill %in% c("blue")]<-"darkslategray1"
symbiontX[4,]$fill[symbiontX[4,]$fill %in% c("blue")]<-"#30C8C4"
symbiontX[5,]$fill[symbiontX[5,]$fill %in% c("blue")]<-"#E57373"
symbiontX[6,]$fill[symbiontX[6,]$fill %in% c("blue")]<-"pink1"
symbiontX[1,]$fill[symbiontX[1,]$fill %in% c("blue")]<-"lightslateblue"
symbiontX$col[symbiontX$col %in% c("blue")]<-"black"

# names
symbiontX <- dna_seg(symbiontX)
dna_segs <- list(symbiontX, symbiontY, symbiontZ)
name <- c("symbiontX", "symbiontY", "symbiontZ")
names(dna_segs) <- name

mid_pos <- middle(dna_seg(symbiontX))
annot_symbiontX <- annotation(x1=mid_pos, 
                               text = symbiontX$name, 
                               rot = c(0), col = c("black"))

# introduce breaks
xlims <- list(c(0, 1330, 1381, 5500),
              c(0, 1450, 1520, 5600),
              c(0, 1400, 1401, 5500),
              c(0, 1400, 1411, 5500),
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL,
              NULL)

# representation
plot_gene_map(dna_segs = dna_segs,
              gene_type = "arrows",
              dna_seg_scale = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
              scale=TRUE, 
              annotations = annot_symbiontX, 
              annotation_height = 0.5, 
              annotation_cex = 0.8,
              scale_cex = 0.8,
              limit_to_longest_dna_seg = FALSE,
              dna_seg_line = c("black"), 
              dna_seg_label_cex = 0.8,
              dna_seg_label_col = "black",
              tree_width = 2, 
              xlims=xlims)
