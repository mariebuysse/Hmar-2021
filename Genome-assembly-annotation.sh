#for the Italian and Spanish isolates (Anna Maria, Tiago, can you complete this part please?)



#for the Israeli isolate
#mapping 
bwa index ref-genome.fasta
bwa mem -t 4 ref-genome.fasta raw-reads_R1.fq raw-reads_R2.fq | samtools view -Sb > outputmapping.bam
samtools sort -@ 4 -o outputsorted.bam outputmapping.bam
samtools index outputsorted.bam 
samtools view -F 2 -b -o unmapped.bam outputsorted.bam
samtools bam2fq unmapped.bam > unmapped.fastq
fastqc unmapped.fastq ##quality checking
samtools flagstat outputmapping.bam > stat.txt ##quality checking

#assembling (SPAdes)
spades.py --12 unmapped.fastq -o draftgenome
quast.py draftgenome_scaffolds.fasta

#process in order to obtain a final assembling genome
ragoo.py draftgenome_scaffolds.fasta ref-genome.fasta

#quality and completeness analyses
quast.py finalgenome_scaffolds.fasta
## ADD MiComplete command line (to be completed by Marie)

#annotation 
prokka --setupdb ./MidiIricVA_current_v2.faa
prokka <FASTA> --outdir <DIR> --genus Midichloria --strain <NAME STRAIN> --cpus 2 --locustag <Midi+NAME STRAIN> --proteins ./MidiIricVA_current_v2.faa
