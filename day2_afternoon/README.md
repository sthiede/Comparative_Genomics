Day 2 Afternoon
===============
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

High-throughput BLAST and pan-genome analysis
---------------------------------------------

This morning we learned how to perform basic genome annotation and comparison using Prokka and ACT. Now we will up the ante and do some more sophisticated comparative genomics analyses! 
First, we will create custom BLAST databases to identify specific antibiotic resistance genes of interest in a set of genomes. 
Second, we will use the large-scale BLAST-based tool LS-BSR to identify the complete antibiotic resistome in our genomes. 
Third, we will move beyond antibiotic resistance, and look at the complete set of protein coding genes in our input genomes. 
Finally, we will go back to ACT to understand the sorts of genomic rearrangements underlying observed variation in gene content.  

For these exercises we will be looking at four closely related Acinetobacter baumannii strains. However, despite being closely related, these genomes have major differences in gene content, as A. baumannii has a notoriously flexible genome! In fact, in large part due to its genomic flexibility, A. baumannii has transitioned from a harmless environmental contaminant to a pan-resistant super-bug in a matter of a few decades. If you are interested in learning more, check out this nature [review](http://www.nature.com/nrmicro/journal/v5/n12/abs/nrmicro1789.html) or [this](http://www.pnas.org/content/108/33/13758.abstract) paper, I published a few years back analyzing the very same genomes you are working with.

Execute the following command to copy files for this afternoon’s exercises to your scratch directory:

```  

cd /scratch/micro612w18_fluxod/username

or

wd

cp -r /scratch/micro612w18_fluxod/shared/data/day2_after/ ./

```

<!---
Changed 23 feb
##Determine which genomes contain beta-lactamase genes
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)
Before comparing full genomic content, lets start by looking for the presence of particular genes of interest. A. baumannii harbors an arsenal of resistance genes, and it would be interesting to know how particular resistance families vary among our 4 genomes. To accomplish this we will use the antibiotic resistance database ([ARDB](http://ardb.cbcb.umd.edu/)). In particular, we are going to extract a set of genes from ARDB that we are interested in probing our genomes for, and create a custom BLAST database to compare against. 
**Get beta-lactamase genes from [ARDB](http://ardb.cbcb.umd.edu/) database**
>i. Run the custom perl script filter_fasta_file.pl to extract genes annotated as beta-lactamases from the full ARDB fasta file. 
The script takes as input: 
1) the ARDB database (resisGenes.pfasta), 
2) a file containing terms to search the database for (fasta_file_keys) and 
3) an output file to contain the subset of sequences that match the text your searching for (ardb_beta_lactam_genes.pfasta).
```
module load bioperl
cd scratch/micro612w18_fluxod/username/day2_after
perl filter_fasta_file.pl resisGenes.pfasta fasta_file_keys ardb_beta_lactam_genes.pfasta
```
>ii. Build BLAST database from fasta file
Run formatdb on the file of beta-lactamases to create a BLAST database. 
formatdb takes as input: 
1) a fasta file of protein or nucleotide sequences (ardb_beta_lactam_genes.pfasta) and 
2) a flag indicating whether to construct a protein or nucleotide database (in this case protein/ -p T).
```
formatdb -i ardb_beta_lactam_genes.pfasta -p T
```
>iii. BLAST A. baumannii proteins against our custom beta-lactamase database
Run BLAST! 
The input parameters are: 
1) the type of blast to use (-p blastp),
2) query sequences (-i Abau_all.pfasta), 
3) the database to search against (-d ardb_beta_lactam_genes.pfasta),
4) the name of a file to store your results (-o bl_blastp_results), 
5) output format (-m 8), 
6) e-value cutoff (-e 1e-20), 
7) number of database sequences to return (-v 1) and 
8) number of database sequences to show alignment for (-b 1).
```
blastall -p blastp -i Abau_all.pfasta -d ardb_beta_lactam_genes.pfasta -o bl_blastp_results -m 8 -e 1e-20 -v 1 -b 1
```
Use less to look at bl_blastp_results.
```
less bl_blastp_results
```
> Question: Experiment with the –m parameter, which controls different output formats that BLAST can produce. 
>iv. Repeat steps i-iii for a different resistance gene class
Use nano to change fasta_file_keys to contain phrase you’d like to search for (e.g. acetyltransferase, carbapenemase)
Run filter_fasta_file.pl to extract genes annotated with your resistance of interest (ROI) from the full ARDB fasta file
```
perl filter_fasta_file.pl resisGenes.pfasta fasta_file_keys ardb_ROI_genes.pfasta
```
Create new BLAST database with formatdb
```
formatdb -i ardb_ROI_genes.pfasta -p T
```
**BLAST!**
```
blastall -p blastp -i Abau_all.pfasta -d ardb_ROI_genes.pfasta -o bl_blastp_results -m 8 -e 1e-20 -v 1 -b 1
```
## Identification of antibiotic resistance genes with [LS-BSR](https://github.com/jasonsahl/LS-BSR) and the [ARDB](http://ardb.cbcb.umd.edu/) database
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)
Next, instead of looking at resistance classes one at a time, lets look at them all in one shot! To do this we will use [LS-BSR](https://peerj.com/articles/332/), which essentially is just a wrapper for doing the same sort of BLASTing we just did in the previous step. BSR stands for BLAST Score Ratio, which refers to what the output is. In particular, for each query gene LS-BSR returns the ratio between: 1) the BLAST score of best hit in target genome and 2) BLAST score of query gene against itself. So, the output is a query by target genome matrix, where the values are between 0 and 1, and indicate the strength of a given queries BLAST hit in the target genome. 
>i. Create a non-redundant list of resistance genes
There is a lot of redundancy in the ARDB (e.g. lots of closely related genes), which would make the output difficult to sort through. Here, we use usearch to select representatives from the database and create a non-redundant gene set! 
We are running usearch with the following parameters: 
1) the clustering algorithm (-cluster_fast), 
2) the files of sequences you want to cluster (resisGenes.pep), 
3) the minimum sequence identity to be included in an existing cluster (-id 0.8), 
4) an output fasta file with reperesentatives (centroids) of each sequence cluster (-centroids resisGenes_nr.pep) and 
5) an output file describing the results of the clustering (-uc resisGenes.uc).
```
> Make sure you are in day2_after directory
cd scratch/micro612w18_fluxod/username/day2_after
> Load relevant Modules
module load usearch 
module load prodigal
> Run usearch to select representatives from the database and create a non-redundant gene set! 
usearch -cluster_fast resisGenes.pep -id 0.8 -centroids resisGenes_nr.pep -uc resisGenes.uc
```
>ii. Run LS-BSR
Change your directory to day2_after:
```
> Make sure you are in day2_after directory
cd /scratch/micro612w18_fluxod/username/day2_after/
```
Run LS-BSR (it will take a few minutes)! 
The input parameters are: a directory with your genomes (-d Abau_genomes) and a fasta file of query genes (-g resisGenes_nr.pep)
```
python /scratch/micro612w18_fluxod/shared/bin/LS-BSR-master/ls_bsr.py -d Abau_genomes/ -g resisGenes_nr.pep
```
>iii. Download LS-BSR output matrix to your own computer for analysis in R
Use scp to get LS-BSR output onto your laptop
```
> Dont forget to change username in the below command
scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w18_fluxod/username/day2_after/bsr_matrix_values.txt ~/Desktop
```
Fire up RStudio and read the matrix:
```
> Make sure you have copied bsr_matrix_values.txt file to your desktop. If not then give the path where bsr_matrix_values.txt is located.
bsr_mat = read.table('~/Desktop/bsr_matrix_values.txt', sep = "\t", row.names = 1, header = TRUE, quote = "")
```
Use head, str, dim, etc. to explore the matrix you read in
iv. Make a heatmap of all the LS-BSR results
Install and load the R library "heatmap3"
Make a heatmap of the complete LS-BSR matrix. Check out the help file to see what the input parameters do, and behold the plethora of other options to customize your heatmaps!
```
heatmap3(bsr_mat, , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```
>v. Subset LS-BSR data to only include genes present in at least one genome 
From the previous step you should have discerned that full LS-BSR matrix is too large to get a useful visualization, so we need to subset it. 
Lets first subset the matrix to focus only on genes present in at least one of our genomes. 
Values in the LS-BSR matrix are between 0 and 1, and represent the sequence identity to the query gene. 
We will arbitrarily say that if a protein have a BLAST score ratio of less then 0.5, then its absent.
```
bsr_mat_subset = bsr_mat[rowSums(bsr_mat > 0.5) > 0,]
```
Make a heatmap of your subset (much better!)
```
heatmap3(bsr_mat_subset, , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```
>vi. Determine the total number of resistance genes present in each genome
We use colSums to count the number of genes with greater than 50% identity to the query
```
colSums(bsr_mat > 0.5)
```
How does the total number of genes vary by altering the percent identity threshold?
>vii. Determine the total number of bla genes in each genome
Next, we will use grepl to pull out genes of interest
```
bla_bsr_mat = bsr_mat[grepl('beta-lactamase', row.names(bsr_mat)) ,]
```
Print out to screen and make a heatmap to explore
>viii. Subset the full matrix to look at genes that are present in only one genome
Get genes present in only one genome
```
unique_bsr_mat = bsr_mat[rowSums(bsr_mat > .5) == 1,]
```
Print out to screen and make a heatmap to explore
## Perform pan-genome analysis with LS-BSR
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)
As a final BLASTing exercise we will use LS-BSR to explore the pan-genome of our A. baumannii. The pan-genome is just a fancy term for the full complement of genes in a set of genomes. 
The way LS-BSR does this is by: 
1) applying prodigal to identify protein coding genes in input genomes, 
2) applying usearch to create non-redundant set of genes and 
3) BLASTing the set of non-redundant genes against the genomes.
>i. Get pan-genome matrix and transfer annotation
Make sure you are on an interactive node, as this will be even more computationally intensive!
```
iflux
```
Change your directory to day2_after
```
> Make sure to change username with your uniqname
cd /scratch/micro612w18_fluxod/username/day2_after/
```
Run LS-BSR! The –u parameter is just a path to where usearch lives on flux.
If you started a new interactive job since you ran LS-BSR, you will need to re-load the required modules for LS-BSR listed above.
```
cd scratch/micro612w18_fluxod/username/day2_after
python /scratch/micro612w18_fluxod/shared/bin/LS-BSR-master/ls_bsr.py -d Abau_genomes/ -u /sw/med/centos7/usearch/8.1/usearch
```
Run the custom perl script transfer_annotations.pl to add annotations to your BSR matrix. The output of this script will be bsr_matrix_values_annot.txt
```
perl transfer_annotations.pl Abau_ECII_PC.fasta Abau_ECII_PC.NR.annot bsr_matrix_values.txt consensus.fasta
```
>ii. Read matrix into R and create heatmap
Use scp to get LS-BSR output onto your laptop
```
> Make sure to change username with your uniqname
scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w18_fluxod/username/day2_after/bsr_matrix_values_annot.txt ~/Desktop
```
Fire up RStudio and read the matrix in
```
bsr_mat_PG = read.table('~/Desktop/bsr_matrix_values_annot.txt', sep = "\t", row.names = 1, header = TRUE, quote = "")
```
Use head, str, dim, etc. to explore the matrix you read in
Make a heatmap for the full matrix
```
heatmap3(as.matrix(bsr_mat_PG), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```
Make a heatmap for variable genes (present in at least one, but not all of the genomes
```
bsr_mat_PG_subset = bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) > 0 & rowSums(bsr_mat_PG > 0.4) < 4 ,] 
heatmap3(as.matrix(bsr_mat_PG_subset), , scale = "none", distfun = function(x){dist(x, method = "manhattan")}, margin = c(10,10), cexCol = 0.85, cexRow = 0.5)
```
>iii. Which genomes are most closely related based upon shared gene content?
We will use the outer function to determine the number of genes shared by each pair of genomes. 
Here we are arbitrarily deciding that a gene is present if the BSR is greater than 0.4. 
Look at the help page for outer to gain additional insight into how this is working.
```
outer(1:4,1:4, FUN = Vectorize(function(x,y){sum(bsr_mat_PG_subset[,x] > 0.4 & bsr_mat_PG_subset[,y] > 0.4)}))
```
>iv. What is the size of the core genome?
Lets first get an overview of how many genes are present in different numbers of genomes (0, 1, 2, 3 or 4) by plotting a histogram. Here, we combine hist with rowSums to accomplish this.
```
hist(rowSums(bsr_mat_PG > 0.4))
```
Next, lets figure out how big the core genome is (e.g. how many genes are common to all of our genomes)?
```
sum(rowSums(bsr_mat_PG > 0.4) == 4)
```
>v. What is the size of the accessory genome?
Lets use a similar approach to determine the size of the accessory genome (e.g. those genes present in only a subset of our genomes).
```
sum(rowSums(bsr_mat_PG > 0.4) < 4 & rowSums(bsr_mat_PG > 0.4) > 0)
```
>vi. What types of genes are unique to a given genome?
So far we have quantified the core and accessory genome, now lets see if we can get an idea of what types of genes are core vs. accessory. Lets start by looking at those genes present in only a single genome. What do you notice about these genes?
```
row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 1,])
```
vii. What is the number of hypothetical genes in core vs. accessory genome?
Looking at unqiue genes we see that many are annotated as “hypothetical”, indicating that the sequence looks like a gene, but has no detectable homology with a functionally characterized gene. Determine the fraction of “hypothetical” genes in unique vs. core. Why does this make sense?
```
sum(grepl("hypothetical" , row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 1,]))) / sum(rowSums(bsr_mat_PG > 0.4) == 1)
sum(grepl("hypothetical" , row.names(bsr_mat_PG[rowSums(bsr_mat_PG > 0.4) == 4,]))) / sum(rowSums(bsr_mat_PG > 0.4) == 4)
```
-->


Determine which genomes contain beta-lactamase genes
----------------------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

Before comparing full genomic content, lets start by looking for the presence of particular genes of interest. A. baumannii harbors an arsenal of resistance genes, and it would be interesting to know how particular resistance families vary among our 4 genomes. To accomplish this we will use the antibiotic resistance database ([ARDB](http://ardb.cbcb.umd.edu/)) and particularly beta-lactamase genes extracted from ARDB. These extracted genes can be found in file ardb_beta_lactam_genes.pfasta, which we will use to generate a Blast database.

>i. Run makeblastdb on the file of beta-lactamases to create a BLAST database. 

makeblastdb takes as input: 

1) an input fasta file of protein or nucleotide sequences (ardb_beta_lactam_genes.pfasta) and 

2) a flag indicating whether to construct a protein or nucleotide database (in this case protein/ -dbtype prot).

```
#change directory to day2_after
d2a


makeblastdb -in ardb_beta_lactam_genes.pfasta -dbtype prot

```

>iii. BLAST A. baumannii protein sequences against our custom beta-lactamase database. 

Run BLAST! 

The input parameters are: 

1) query sequences (-query Abau_all.pfasta), 

2) the database to search against (-db ardb_beta_lactam_genes.pfasta), 

3) the name of a file to store your results (-out bl_blastp_results), 

4) output format (-outfmt 6), 

5) e-value cutoff (-evalue 1e-20), 

6) number of database sequences to return (-max_target_seqs 1)


```

blastp -query Abau_all.pfasta -db ardb_beta_lactam_genes.pfasta -out bl_blastp_results -outfmt 6 -evalue 1e-20 -max_target_seqs 1

```

Use less to look at bl_blastp_results.

```
less bl_blastp_results
```

> Question: Experiment with the –outfmt parameter, which controls different output formats that BLAST can produce. 

Identification of antibiotic resistance genes with [ARIBA](https://github.com/sanger-pathogens/ariba) directly from paired end reads
----------------------------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

ARIBA, Antimicrobial Resistance Identification By Assembly is a tool that identifies antibiotic resistance genes by running local assemblies. The input is a FASTA file of reference sequences (can be a mix of genes and noncoding sequences) and paired sequencing reads. ARIBA reports which of the reference sequences were found, plus detailed information on the quality of the assemblies and any variants between the sequencing reads and the reference sequences.

ARIBA is compatible with various databases and also contains an utility to download different databases such as: argannot, card, megares, plasmidfinder, resfinder, srst2_argannot, vfdb_core. Today, we will be working with the [card](https://card.mcmaster.ca/) database, which has been downloaded and placed in /scratch/micro612w18_fluxod/shared/out.card.prepareref/ directory.

<!---
Note: There is an issue with downloading the database. They are in a process to fix the broken CARD database link issue. For now, I am using my own downloaded database.
>i. Run this command to download CARD database
```
/nfs/esnitkin/bin_group/anaconda3/bin/ariba getref card out.card
```
>ii. Prepare this downloaded card database for ARIBA
```
/nfs/esnitkin/bin_group/anaconda3/bin/ariba prepareref -f out.card.fa -m out.card.tsv out.card.prepareref
```
-->

>i. Run ARIBA on input paired-end fastq reads for resistance gene identification. 

The fastq reads are placed in Abau_genomes_fastq directory. Enter interactive flux session, change directory to day2_after workshop directory and run the below four commands to start ARIBA jobs in background.

```
iflux

cd /scratch/micro612w18_fluxod/username/day2_after

or 

d2a

#Load dependency

module load cd-hit

#ARIBA commands

/nfs/esnitkin/bin_group/anaconda3/bin/ariba run --force /scratch/micro612w18_fluxod/shared/out.card.prepareref/ Abau_genomes_fastq/AbauA_genome.1.fastq.gz Abau_genomes_fastq/AbauA_genome.2.fastq.gz AbauA_genome &

/nfs/esnitkin/bin_group/anaconda3/bin/ariba run --force /scratch/micro612w18_fluxod/shared/out.card.prepareref/ Abau_genomes_fastq/AbauB_genome.1.fastq.gz Abau_genomes_fastq/AbauB_genome.2.fastq.gz AbauB_genome &

/nfs/esnitkin/bin_group/anaconda3/bin/ariba run --force /scratch/micro612w18_fluxod/shared/out.card.prepareref/ Abau_genomes_fastq/AbauC_genome.1.fastq.gz Abau_genomes_fastq/AbauC_genome.2.fastq.gz AbauC_genome &

/nfs/esnitkin/bin_group/anaconda3/bin/ariba run --force /scratch/micro612w18_fluxod/shared/out.card.prepareref/ Abau_genomes_fastq/ACICU_genome.1.fastq.gz Abau_genomes_fastq/ACICU_genome.2.fastq.gz ACICU_genome &

```

The "&" in the above commands(at the end) is a little unix trick to run commands in background. You can run multiple commands in background and make full use of parallel processing. You can check the status of these background jobs by typing:

```
jobs
```

>ii. Run ARIBA summary function to generate a report summary.

ARIBA has a summary function that summarises the results from one or more sample runs of ARIBA and generates an output report with various level of information determined by -preset parameter. The parameter "-preset minimal" will generate a minimal report showing only the presence/absence of resistance genes whereas "-preset all" will output all the extra information related to each database hit such as reads and reference sequence coverage, variants and their associated annotations(if the variant confers resistance to an Antibiotic) etc.

```

/nfs/esnitkin/bin_group/anaconda3/bin/ariba summary --preset minimal Abau_genomes_ariba_minimal_results *_genome/*.tsv

/nfs/esnitkin/bin_group/anaconda3/bin/ariba summary --preset all Abau_genomes_ariba_all_results *_genome/*.tsv

```

ARIBA summary generates three output:

1. Abau_genomes_ariba*.csv file that can be viewed in your favourite spreadsheet program.
2. Abau_genomes_ariba*.phandango.{csv,tre} that allow you to view the results in [Phandango](http://jameshadfield.github.io/phandango/#/). They can be drag-and-dropped straight into Phandango.

Download these files with Cyberduck or the scp command. 

Perform pan-genome analysis with [Roary](https://sanger-pathogens.github.io/Roary/)
----------------------------------------

Roary is a pan genome pipeline, which takes annotated assemblies in GFF3 format and calculates the pan genome. The pan-genome is just a fancy term for the full complement of genes in a set of genomes. 

The way Roary does this is by: 
1) Roary gets all the coding sequences from GFF files, convert them into protein, and create pre-clusters of all the genes, 
2) Then, using BLASTP and MCL, Roary will create gene clusters, and check for paralogs. and 
3) Finally, Roary will take every isolate and order them by presence/absence of genes.

>i. Generate pan-genome matrix using Roary and GFF files

Make sure you are on an interactive node, as this will be even more computationally intensive!

```
iflux
```

Change your directory to day2_after

```

> Make sure to change username with your uniqname

cd /scratch/micro612w18_fluxod/username/day2_after/

or 

d2a

```

Load all the required dependencies and run roary on GFF files placed in Abau_genomes_gff folder.

```
module load samtools
module load bedtools2
module load cd-hit
module load ncbi-blast
module load mcl
module load parallel
module load mafft
module load fasttree
module load perl-modules
module load R
module load roary

#Run roary
roary -p 4 -f Abau_genomes_roary_output -r -n -v Abau_genomes_gff/*.gff 
```

The above roary command will run pan-genome pipeline on gff files placed in Abau_genomes_gff(-v) using 4 threads(-p), save the results in an output directory Abau_genomes_roary_output(-f), generate R plots using .Rtab output files and align core genes(-n)

Change directory to Abau_genomes_roary_output to explore the results.

```
cd Abau_genomes_roary_output

ls
```

Output files:

1. summary_statistics.txt: This file is an overview of your pan genome analysis showing the number of core genes(present in all isolates) and accessory genes(genes absent from one or more isolates or unique to a given isolate). 

2. gene_presence_absence.csv: This file contain detailed information about each gene including their annotations which can be opened in any spreadsheet software to manually explore the results. It contains plethora of information such as gene name and their functional annotation, whether a gene is present in a genome or not, minimum/maximum/Average sequence length etc.

3. gene_presence_absence.Rtab: This file is similar to the gene_presence_absence.csv file, however it just contains a simple tab delimited binary matrix with the presence and absence of each gene in each sample. It can be easily loaded into R using the read.table function for further analysis and plotting. The first row is the header containing the name of each sample, and the first column contains the gene name. A 1 indicates the gene is present in the sample, a 0 indicates it is absent.

4. core_gene_alignment.aln: a multi-FASTA alignment of all of the core genes that can be used to generate a phylogenetic tree.

<!--
#Plots are not very useful. Seems like a waste of time.
>ii. Generate a phylogenetic tree and plot pan-genome matrix.
we will use core_gene_alignment.aln multi-fasta core gene alignment as an input to generate a phylogenetic tree using FastTree tool. 
This tree along with the pan-genome matrix can then be used to plot some nice plots. Roary comes with a python script called roary_plots.py that can be used for visualizing pan-genome analysis results. 
```
module load fasttree
FastTree core_gene_alignment.aln > core_gene_alignment.tree
module load python-anaconda3/latest
python /scratch/micro612w18_fluxod/shared/bin/roary/contrib/roary_plots/roary_plots.py core_gene_alignment.tree gene_presence_absence.csv
```
-->

>ii. Explore pan-genome matrix gene_presence_absence.csv and gene_presence_absence.Rtab using R

Note:plots generated by roary_plots.py doesn't seem to be very useful and is completely a waste of time. query_pan_genome script provided by roary doesn't work and generates an empty result which seems like a bug and the link to the issue raised for this bug can be found [here](https://github.com/sanger-pathogens/Roary/issues/298) 



Perform genome comparisons with [ACT](http://www.sanger.ac.uk/science/tools/artemis-comparison-tool-act)
-------------------------------------
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day2_afternoon/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

In the previous exercises we were focusing on gene content, but losing the context of the structural variation underlying gene content variation (e.g. large insertions and deletions). 
Here we will use ACT to compare two of our genomes (note that you can use ACT to compare more than two genomes if desired). 

i. Create ACT alignment file with BLAST

As we saw this morning, to compare genomes in ACT we need to use BLAST to create the alignments. We will do this on flux.

```

cd scratch/micro612w18_fluxod/username/day2_after
blastall -p blastn -i ./Abau_genomes/AbauA_genome.fasta -d ./Abau_BLAST_DB/ACICU_genome.fasta -m 8 -e 1e-20 -o AbauA_vs_ACICU.blast

```

>ii. Read in genomes, alignments and annotation files

Use scp or cyberduck to transfer Abau_ACT_files folder onto your laptop


1. Abau_genomes/AbauA_genome.fasta 
2. Abau_genomes/ACICU_genome.fasta 
3. AbauA_vs_ACICU.blast 
4. Abau_ACT_files/AbauA_genome_gene.gff 
5. Abau_ACT_files/ACICU_genome_gene.gff


>iii. Explore genome comparison and features of ACT

Read in genomes and alignment into ACT

```

Go to File -> open 
Sequence file 1  = ACICU_genome.fasta 
Comparison file 1  = AbauA_vs_ACICU.blast
Sequence file 2  = AbauA_genome.fasta

```

Before we use annotation present in genbank files. Here we will use ACT specific annotation files so we get some prettier display (resistance genes = red, transposable elements = bright green)  

```

Go to File -> ACICU_genome.fasta -> Read an entry file = ACICU_genome_gene.gff

Go to File -> AbauA_genome.fasta -> Read an entry file = AbauA_genome_gene.gff

```

Play around in ACT to gain some insight into the sorts of genes present in large insertion/deletion regions. 
See if you can find: 
1) differences in phage content, 
2) membrane biosynthetic gene cluster variation and 
3) antibiotic resistance island variation.

