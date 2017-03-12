# Day 1 Morning
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

## Getting your data onto Flux and setting up environment variable

**Log in to Flux**

```
ssh user@flux-login.arc-ts.umich.edu
```

<!-- **Set up your .bashrc file so your environment is all set for genomic analysis!** -->

**Setting up environment variables in .bashrc file so your environment is all set for genomic analysis!**

Environment variables are the variables/values that describe the environment in which programs run in. All the programs and scripts on your unix system use these variables for extracting information such as: What is my current working directory?, Where are temporary files stored?, Where are perl/python libraries?, Where is Blast installed? etc. 

In addition to environment variables that are set up by system administators, each user can set their own environment variables to customize their experience. This may sound like something super advanced that isn't relevent to beginners, but that's not true! Some examples of ways that we will use environment variables in the class are: 1) create shortcuts for directories that you frequently go to, 2) tell unix where frequently used programs live, so you don't have to put the full path name each time you use it and 3) setup a shortcut for getting on a cluster node, so that you don't have to write out the full command each time.

One way to set your environment variables would be to manually set up these variables everytime you log in, but this would be extremely tedious and inefficient. So, Unix has setup a way around this, which is to put your environment variable assignments in special files called .bashrc or .bash_profile. Every user has one or both of these files in thier home directory, and what's special about them is that the commands in them are executed every time you login. So, if you simply set your environmental variable assignments in one of these files, your environment will be setup just the way you want it each time you login!

All the softwares/tools that we need in this workshop are installed in a directory "/scratch/micro612w17_fluxod/shared/bin/" and we want the shell to look for these installed tools in this directory. For this, We will save the full path to these tools in an environment variable PATH.

>i. Make a backup copy of bashrc file in case something goes wrong. 
	
```

cp ~/.bashrc ~/bashrc_backup

Note: "~/" represents your home directory. On flux, these means /home/username

```
	
>ii. Open ~/.bashrc file using any text editor and add the following lines to your .bashrc file.

<details>
  <summary>Click to expand entries</summary>
  
```
## Micro612 Workshop ENV

#Aliases
alias iflux='qsub -I -V -l nodes=1:ppn=1,mem=4000mb,walltime=1:00:00:00 -q fluxod -l qos=flux -A micro612w17_fluxod'

# Flux Modules
module load python-anaconda2/latest
module load perl-modules

# Perl Libraries
export PERL5LIB=/scratch/micro612w17_fluxod/shared/bin/PAGIT/lib:/scratch/micro612w17_fluxod/shared/bin/vcftools_0.1.12b/perl:$PERL5LIB
export PERL5LIB=/scratch/micro612w17_fluxod/shared/perl_libs:$PERL5LIB

# Bioinformatics Tools
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/mauve_snapshot_2015-02-13/linux-x64/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/blast/bin/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/vcftools_0.1.12b/perl/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/tabix-0.2.6/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/bwa-0.7.12/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/Trimmomatic/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/bcftools-1.2/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/samtools-1.2/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/sratoolkit/bin/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/Spades/bin/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/FastQC/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/GenomeAnalysisTK-3.3-0/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/picard-tools-1.130/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/qualimap_v2.1/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/vcftools_0.1.12b/bin/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/snpEff/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/PAGIT/ABACAS/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/blast-2.2.26/bin/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/quast/
export PATH=$PATH:/scratch/micro612w17_fluxod/shared/bin/MUMmer3.23/

```

</details>

The above environment settings will set a shortcut "iflux" for entering interactive flux session, call necessary flux modules and perl libraries required by certain tools and finally sets the path for bioinformatics programs that we will run during the workshop.

>iii. Save the file and Source .bashrc file to make these changes permanent.

```

source ~/.bashrc

```

>iv. Check if the $PATH environment variable is updated

```

echo $PATH

```


<!-- Check the dependencies Pending
tree file system Pending
-->


**Power of Unix


In software carpentry, you learned working with shell and automating simple tasks using basic unix commands. Lets see how some of these commands can be employed in genomics analysis while exploring various file formats that we use in day to day analysis. For these session, we will be try to explore three different type of bioinformatics file formats: 

fasta: used for representing either nucleotide or peptide sequences

gff: used for describing genes and other features of DNA, RNA and protein sequences

fastq: used for storing biological sequence / sequencing reads (usually nucleotide sequence) and its corresponding quality scores

> Execute the following commands to copy files for this morning’s exercises to your home directory: 

```

cd /scratch/micro612w17_fluxod/username
cp -r /scratch/micro612w17_fluxod/shared/data/day1_morn/ ./
cd /scratch/micro612w17_fluxod/username/day1_morn/
ls

```

> Question: In the homework assignment, you downloaded genome assembly fasta files and ran a shell script to count the contigs. Lets say you want to find out the combined length of genome in each of these files. This can be achieved by running a simple unix command comprised of grep, sed and awk. The key here is knowing the features of fasta file format such as: each sequence is preceded by a fasta header that starts with ">", types of bases that a nucleotide sequence represents (A,T,G,C,N) and each line is seperated by a new line character ("\n"). To achieve this, we will use grep to match only those lines that doesn't start with ">" (remember grep -v option to ignore lines), use sed to remove characters that match "N" or "n" which represents unknown bases and finally use awk to count the remaining characters. We can use unix pipe "|" to pass the output of one command to another for further processing.


<details>
  <summary>Solution</summary>
```
grep -v '^>' filename.fasta | sed 's/[N,n]//g' | awk -F '\n' '{sum += length} END {print sum}'

Note:
- The above command will ignore fasta header lines (starts with ">"), replace characters "N" and "n" and uses awk to count the remaining characters. 
- The sign "^" inside the grep pattern represents any pattern that starts with "^>" and -v asks grep to ignore those lines.
- Use "|" to pass these lines to sed. sed stands for stream editor which parses, transforms and replaces text. Here, we are removing the characters "N" or "n" and keeping only "A,T,G,C" bases
- awk consists of three blocks: The first block (-F '\n') tells awk how each line is seperated from each other using a field seperator, the second block tells awk to keep counting characters in a line (using awk's default option "length") and save the count in a variable "sum" and when it runs through all the lines in a stream, the third block prints the value of sum which represents total bases in a fasta file.

```

</details>

## Quality Control using [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ "FastQC homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day1_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)



As soon as you receive your sample data from sequencing centre, the first thing you do is check its quality using a quality control tool such as FastQC. But before carrying out extensive QC, you can run a bash "one-liner" to get some basic statistics about the raw reads.These one-liners are great examples for how a set of simple (relatively) Unix commandscan be piped together to do really useful things.

Run the following command to print total number of reads in each file, total number of unique reads, percentage of unique reads, most abundant sequence(useful to find adapter sequences or contamination), its frequency, and frequency of that sequence as a proportion of the total reads.

```
for i in *.gz; do zcat $i | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(!max||count[read]>max) {max=count[read];maxRead=read};if(count[read]==1){unique++}};print total,unique,unique*100/total,maxRead,count[maxRead],count[maxRead]*100/total}'; done
```

You can find more of such super useful bash one-liners at Stephen Turner's github [page.](https://github.com/stephenturner/oneliners)

Now we will run FastQC on these raw data to assess its quality. FastQC is a quality control tool that reads in sequence data in a variety of formats(fastq, bam, sam) and can either provide an interactive application to review the results or create an HTML based report which can be integrated into any pipeline. It is generally the first step that you take upon receiving the sequence data from sequencing facility to get a quick sense of its quality and whether it exhibits any unusual properties (e.g. contamination or unexpected biological features)

>ii. In your day1_morn directory, create a new directory for saving FastQC results.

```
mkdir Rush_KPC_266_FastQC_results
mkdir Rush_KPC_266_FastQC_results/before_trimmomatic
```

>iii. Verify that FastQC is in your path by invoking it from command line.

```
fastqc -h
```

FastQC can be run in two modes: "command line" or as a GUI (graphical user interface). We will be using command line version of it.

>iv. Get an interactive cluster node to start running programs

```
qsub -I -V -l nodes=1:ppn=1,mem=4000mb,walltime=1:00:00:00 -q fluxod -l qos=flux -A micro612w17_fluxod
```

>v. Run FastQC to generate quality report of sequence reads.

```
fastqc -o Rush_KPC_266_FastQC_results/before_trimmomatic/ Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz --extract
```

This will generate two results directory, Rush_KPC_266_1_combine_fastqc and Rush_KPC_266_2_combine_fastqc in output folder provided with -o flag. 
The summary.txt file in these directories indicates if the data passed different quality control tests in text format.
You can visualize and assess the quality of data by opening html report in a local browser.

>vi. Exit your cluster node so you don’t waste cluster resources and $$$!

>vii. Download the FastQC report to your home computer to examine

```
sftp username@flux-login.arc-ts.umich.edu
cd /scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/before_trimmomatic/
get Rush_KPC_266_1_combine_fastqc.html
get Rush_KPC_266_2_combine_fastqc.html

or use scp.

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/before_trimmomatic/*.html /path-to-local-directory/
```

The analysis in FastQC is broken down into a series of analysis modules. The left hand side of the main interactive display or the top of the HTML report show a summary of the modules which were run, and a quick evaluation of whether the results of the module seem entirely normal (green tick), slightly abnormal (orange triangle) or very unusual (red cross). 

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day1_morning/1.png)

Lets first look at the quality drop(per base sequence quality graph) at the end of "Per Base Sequence Quality" graph. This degredation of quality towards the end of reads is commonly observed in illumina samples. The reason for this drop is that as the number of sequencing cycles performed increases, the average quality of the base calls, as reported by the Phred Scores produced by the sequencer falls. 

Next, lets check the overrepresented sequences graph and the kind of adapters that were used for sequencing these samples (Truseq or Nextera) which comes in handy while indicating the adapter database during downstream filtering step.

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day1_morning/2.png)

Check out [this](https://sequencing.qcfail.com/articles/loss-of-base-call-accuracy-with-increasing-sequencing-cycles/) for more detailed explaination as to why quality drops with increasing sequencing cycles.

> [A video FastQC walkthrough created by FastQC developers](https://www.youtube.com/watch?v=bz93ReOv87Y "FastQC video") 

## Quality Trimming using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic "Trimmomatic Homepage")
[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day1_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)

Filtering out problematic sequences within a dataset is inherently a trade off between sensitivity (ensuring all contaminant sequences are removed) and specificity (leaving all non-contaminant sequence data intact). Adapter and other technical contaminants can potentially occur in any location within the reads.(start, end, read-through, partial adapter sequences)

Trimmomatic is a tool that tries to search these potential contaminant/adapter sequence within the read at all the possible locations. It takes advantage of the added evidence available in paired-end dataset. In paired-end data, read-through/adapters can occur on both the forward and reverse reads of a particular fragment in the same position. Since the fragment was entirely sequenced from both ends, the non-adapter portion of the forward and reverse reads will be reverse-complements of each other. This strategy of searching for contaminant in both the reads is called 'palindrome' mode. 
For more information on how Trimmomatic tries to achieve this, Please refer [this](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) manual.

Now we will run Trimmomatic on these raw data to remove low quality reads as well as adapters. 

>i. Get an interactive cluster node to start running programs

```
qsub -I -V -l nodes=1:ppn=4,mem=16000mb,walltime=1:00:00:00 -q fluxod -l qos=flux -A micro612w17_fluxod
```

Change your directory to day1_morn

```
cd /scratch/micro612w17_fluxod/username/day1_morn/
```

>ii. Create these output directories in your day1_morn folder to save trimmomatic results

```
mkdir Rush_KPC_266_trimmomatic_results
```

>iii. Load latest version of java and try to invoke trimmomatic from command line.

```
module load lsa java/1.8.0

java -jar /scratch/micro612w17_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar –h
```

>iv. Run the below trimmomatic commands on raw reads.

```
time java -jar /scratch/micro612w17_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/forward_unpaired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w17_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:15 MINLEN:40 HEADCROP:0
```

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day1_morning/trimm_parameters.png)

First, Trimmomatic searches for any matches between the reads and adapters. Short sections(2 bp) of each adapters determined by seed misMatch parameter are tested in each possible position within the reads. If it finds a perfect match, the full alignment is scored to see if the adaptor is truly present. The advantage here is that the full alignment is calculated only when there is a perfect seed match which results in considerable efficiency gains. So, When it finds a match, it moves forward with full alignment and when the match reaches 10 bp determined by simpleClipThreshold, it finally trims off the adapter from reads.  

Quoting Trimmomatic:

"'Palindrome' trimming is specifically designed for the case of 'reading through' a short fragment into the adapter sequence on the other end. In this approach, the appropriate adapter sequences are 'in silico ligated' onto the start of the reads, and the combined adapter+read sequences, forward and reverse are aligned. If they align in a manner which indicates 'read- through' i.e atleast 30 bp match, the forward read is clipped and the reverse read dropped (since it contains no new data)."

>v. Now create new directories in day1_morn folder and Run FastQC on these trimmomatic results.

```
mkdir Rush_KPC_266_FastQC_results/after_trimmomatic

fastqc -o Rush_KPC_266_FastQC_results/after_trimmomatic/ Rush_KPC_266_trimmomatic_results/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results/reverse_paired.fq.gz --extract
```

Get these html reports to local system.
```
sftp username@flux-login.arc-ts.umich.edu
cd /scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic/
get forward_paired.fq_fastqc.html
get reverse_paired.fq_fastqc.html

or use scp 

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic/*.html /path-to-local-directory/
```

![alt tag](https://github.com/alipirani88/Comparative_Genomics/blob/master/_img/day1_morning/3.png)

After running Trimmomatic, you should notice that the sequence quality improved and now doesn't contain any contaminants/adapters.

Next, take a look at the per base sequence content graph, and notice that the head bases(~9 bp) are slightly imbalanced. In a perfect scenario, each nucleotide content should run parallel to each other, and should be reflective of the overall A/C/T/G content of your input sequence. 

Quoting FastQC:
	"It's worth noting that some types of library will always produce biased sequence composition, normally at the start of the read. Libraries produced by priming using random hexamers (including nearly all RNA-Seq libraries) and those which were fragmented using transposases inherit an intrinsic bias in the positions at which reads start. This bias does not concern an absolute sequence, but instead provides enrichment of a number of different K-mers at the 5' end of the reads. Whilst this is a true technical bias, it isn't something which can be corrected by trimming and in most cases doesn't seem to adversely affect the downstream analysis. It will however produce a warning or error in this module."

This doesn't look very bad but you can remove the red cross sign by trimming these imbalanced head bases using HEADCROP:9 flag in the above command.

>vi. Lets Run trimmomatic again with headcrop 9 and save it in a different directory called Rush_KPC_266_trimmomatic_results_with_headcrop/

```
mkdir Rush_KPC_266_trimmomatic_results_with_headcrop/

time java -jar /scratch/micro612w17_fluxod/shared/bin/Trimmomatic/trimmomatic-0.33.jar PE Rush_KPC_266_1_combine.fastq.gz Rush_KPC_266_2_combine.fastq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/forward_unpaired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_unpaired.fq.gz ILLUMINACLIP:/scratch/micro612w17_fluxod/shared/bin/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10:8:true SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:9
```

>vii. Run FastQC 'one last time' on updated trimmomatic results with headcrop and check report on your local computer

```
mkdir Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/
fastqc -o Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/ --extract -f fastq Rush_KPC_266_trimmomatic_results_with_headcrop/forward_paired.fq.gz Rush_KPC_266_trimmomatic_results_with_headcrop/reverse_paired.fq.gz
```
Download the reports again and see the difference.
```
sftp username@flux-login.arc-ts.umich.edu
cd /scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/
get forward_paired.fq_fastqc.html
get reverse_paired.fq_fastqc.html

or use scp

scp username@flux-xfer.arc-ts.umich.edu:/scratch/micro612w17_fluxod/username/day1_morn/Rush_KPC_266_FastQC_results/after_trimmomatic_headcrop/*.html /path-to-local-directory/
```

The red cross sign disappeared!

Lets have a look at one of the Bad Illumina data example [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

[[back to top]](https://github.com/alipirani88/Comparative_Genomics/blob/master/day1_morning/README.md)
[[HOME]](https://github.com/alipirani88/Comparative_Genomics/blob/master/README.md)
