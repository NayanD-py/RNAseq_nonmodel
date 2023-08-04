# Tutorial for de novo transcriptome assembly, annotation and differential gene expression analysis in a non-model organism. 
  
This repository is a usable, publicly available tutorial for analyzing differential expression data. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use `sbatch scriptname` after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bioinformatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.  
  
### Contents  
1. [Introduction](#1-introduction)  
2. [Quality Control](#2-quality-control)   
3. [Assembling Transcriptomes](#3-assembling-transcriptomes)  
4. [Identifying the Coding Regions](#4-identifying-the-coding-regions)  
5. [Determining and Removing Redundant Transcripts](#5-determining-and-removing-redundant-transcripts)
6. [Evaluating the Assembly](#6-evaluating-the-assembly)
7. [Functional Annotation](#7-functional-annotation) 
8. [Quantifying Gene Expression](#8-quantifying-gene-expression)
9. [Differential Expression Analysis in R](#9-Differential-Expression-Analysis-In-R)


## 1. Introduction  
  
The goal of this tutorial is to guide you through a differential gene expression analysis using RNA-seq data in a **non-model organism** for which there is no reference genome available. When a reference genome (or transcriptome) and accompanying functional annotation is available, it is usually preferable to map reads against that resource to do differential expression and downstream functional analysis. 

This isn't always the case, however, so this tutorial includes steps for assembling and annotating a reference transcriptome from RNA-seq data gathered for a differential gene expression analysis, with a resulting workflow that looks like this: 

![ out line ](/images/outline_wide.png)  

The data we'll use here is from a study of needle abscission in the Eastern larch, a deciduous conifer and decidedly **non-model organism**. The study examines gene expression in abscission zone tissue (between the needle and the stem) during the process of needle loss at three time points in autumn, but for this tutorial we will only compare two of them. The leaf samples we'll analyze here were collected from three trees at two time points each (roughly one month apart). This results in six sequence libraries. The libraries were paired-end sequenced with 76bp reads. 

  Sample   |   Location   |   Time point  |   Replicate tree
  ----  |   ----   |   ----  |   ----
  K21   |  Killingworth  |   2  |   1 
  K22   |  Killingworth  |   2  |   2   
  K23   |  Killingworth  |   2  |   3
  K31   |  Killingworth  |   3  |   1 
  K32   |  Killingworth  |   3  |   2
  K33   |  Killingworth  |   3  |   3
  

### Cloning the workflow 

In this workflow we have separated each step into folders. In each folder is one or more scripts for each step in the analysis. Results are written to the same directory. 

To work through the tutorial first clone the repository: 
```bash
git clone https://github.com/CBC-UCONN/RNAseq_nonmodel.git
```

Once you clone the repository you can see the following folder structure:  

```  
Eastern_larch/
├── 01_Raw_Reads
├── 02_Quality_Control
├── 03_Assembly
├── 04_Coding_Regions
├── 05_Clustering
├── 06_RNAQuast
├── 07_EnTAP
├── 08_Counts
├── 09_R_analysis
```  
   
###  SLURM scripts   
This tutorial is tailored to the University of Connecticut's Xanadu computing cluster. Before beginning the tutorial, we'll need to understand a few aspects of the cluster. When first logging into Xanadu from your local terminal, you will be connected to a __submit node__. The submit node provides an interface for interacting with the cluster. Heavy computational work is not done on submit nodes, but rather on __compute nodes__. Never, under any circumstance do substantial computation on a submit node. You may inspect and move files around, or do light editing of scripts, but if you are running a command that uses more than a single processor, or takes more than a couple minutes you should be using a compute node. If you do use a submit node, your process could be killed, all of your work lost, and you may receive a grumpy e-mail or slack message from a cluster administrator. 

Access to compute nodes is managed by [SLURM](https://slurm.schedmd.com/documentation.html), software that allocates computational resources to users of the cluster. Jobs we will run in this tutorial will be processed in __batch mode__, meaning that we pass a script to SLURM, which sends it to a compute node to be run (in contrast to __interactive mode__, where we ask for computational resources and then use them to manually run code). Batch scripts each need a header section that tells SLURM what resources the job needs. Header sections look like this:  

```bash
#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
```  

Following the header section is the code that is to be run. A batch job can be submitted by the command `sbatch myscript.sh`. This tutorial is not intended to teach you the details of the Xanadu cluster and SLURM, so before analyzing your own data, it would be useful to be familiar with the topics covered in the [Xanadu tutorial](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/xanadu/) and our [guidance on resource allocation requests in SLURM](https://github.com/CBC-UCONN/CBC_Docs/wiki/Requesting-resource-allocations-in-SLURM).  


### Obtaining the data 

These data files are only available through Xanadu cluster as they are unpublished. The scripts assume the data are located in the `01_Raw_Reads/` directory. To avoid redundancy on the cluster, rather than copying the data there, you can create __symlinks__ to the read files. Symlinks are "symbolic links"--that is, they are file handles that point at the original files rather than containing the original data. 

In the `01_Raw_Reads/` directory we have written a script ([raw_data_symlinks.sh](/01_Raw_Reads/raw_data_symlinks.sh)) to create these symlinks. You can run this script from the `01_Raw_Reads` directory using `sbatch raw_data_symlinks.sh`.  

Although the script uses a loop, the basic syntax for creating a symlink is `ln -s target_file symlink_name`
  
After running the script, the `01_Raw_Reads` folder will look something like this:  

```
01_Raw_Reads/
├── K21_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K21_R1.fastq.gz
├── K21_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K21_R2.fastq.gz
├── K22_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K22_R1.fastq.gz
├── K22_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K22_R2.fastq.gz
├── K23_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K23_R1.fastq.gz
├── K23_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K23_R2.fastq.gz
├── K31_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K31_R1.fastq.gz
├── K31_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K31_R2.fastq.gz
├── K32_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K32_R1.fastq.gz
├── K32_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K32_R2.fastq.gz
├── K33_R1.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K33_R1.fastq.gz
├── K33_R2.fastq.gz -> /UCHC/PublicShare/CBC_Tutorials/rnaseq_nonmodel_eastern_larch/raw_data/K33_R2.fastq.gz
├── raw_data_symlinks_1979294.err
├── raw_data_symlinks_1979294.out
└── raw_data_symlinks.sh
```

The "->" indicates that the filenames are symlinks, and gives the location of the linked file. 

### Familiarizing yourself with the raw read data

The reads with which we will be working have been sequenced using the [Illumina] platform (https://www.illumina.com/techniques/sequencing.html). We'll assume here that you are familiar with the sequencing technology. Let's have a look at the content of one of our reads, which are in the `fastq` format. They are additionally compressed using `gzip`, so we'll use `zcat`, a pipe (`|`) and `head` to decompress and view the first sequence record:

```bash
zcat K32_R1.fastq | head -n 4
```

`zcat` decompresses and writes the fastq file in plain text. The pipe (`|`) captures the output and sends it to `head -n 4` which will print the first four lines to the screen and then stops the process: 

```
@NS500402:381:HH3NFBGX9:1:11101:2166:1038 1:N:0:CGCTCATT+AGGCTATA
AGAACTCGAAACTAAACGTGGACGTGNTNNTATAAACNNANACNAATCCATCGCCGGTTNNCNTATNNNNNNNNNN
+
AAAAAEEEEEEEEEEEEEEEEEEEEE#E##EEEEEEE##E#EE#EEEE6EEEEEEEEEE##A#EAE##########
```

The first line is the read ID which always begins with "@". The second gives the nucleotide sequence of the read. The third is a comment line, which always begins with "+" and is usually otherwise empty. The fourth line gives the [phred-scaled quality score](https://en.wikipedia.org/wiki/Phred_quality_score) for each base call in the read. The scores are [encoded using ASCII characters](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm) and give the estimated probability that the called base is in error. 
   
   
## 2. Quality Control

### Quality control of Illumina reads using FastQC/MultiQC and Trimmommatic

The first step in analyzing the sequence data is to evaluate its quality. The first thing we'll do is run a pair of programs, `FastQC` and `MultiQC` to assess some basic aspects of the quality of the data that came off the sequencer. `FastQC` generates reports for each fastq file and `MultiQC` aggregates the individual reports into a single HTML file to make it easy to see look them all over quickly. Then we'll run `Trimmomatic` to trim low quality sequence and adapter contamination from our reads. Finally we'll run  FastQC/MultiQC again to see how our data have improved. 

To run FastQC/MultiQC, move to the directory `02_Quality_Control` and run the script `qc_raw.sh` by entering `sbatch qc_raw.sh` on the command line. 

A `FastQC` command line simply invokes the program, the output directory you want to use, and any number of fastq files, like this:

```bash
fastqc --outdir ./raw_fastqc/ ../01_Raw_Reads/K22_R1.fastq.gz ../01_Raw_Reads/K22_R2.fastq.gz
```

`MultiQC` is run like this:

```bash
multiqc --outdir raw_multiqc ./raw_fastqc/
```

This job will take about 30 minutes to complete. When it has completed you should download the HMTL report `multiqc_report.html` found in the newly created `raw_multiqc` directory. You can do this using `scp` or a GUI file transfer program such as Cyberduck or FileZilla. See the [CBC guidance on file transfers](https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/data-transfer-2/) for more details. 

Overall our data looks pretty good. FastQC flags some statistics as "Fail" but that's because the expectations are calibrated to whole-genome sequencing. In RNA-seq, some non-random base composition at the beginning of the sequences, as well as some overrepresented sequences are expected. For detailed explanation of each module in `FastQC` see [the documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/). 

Although our data looks good, we'll still run `Trimmomatic` to remove low quality and residual adapter sequence. This is critical for transcriptome assembly. You can run the `Trimmomatic` script from the `02_Quality_Control` directory by entering `sbatch trimmomatic.sh` on the command line. This script will take about 2 hours to run. To run Trimmomatic on a single sample looks like this:

```bash
SAM=K21
java -jar $Trimmomatic PE -threads 4 \
        ../01_Raw_Reads/${SAM}_R1.fastq.gz \
        ../01_Raw_Reads/${SAM}_R2.fastq.gz \
        trim_${SAM}_R1.fastq.gz singles_trim_${SAM}_R1.fastq.gz \
        trim_${SAM}_R2.fastq.gz singles_trim_${SAM}_R2.fastq.gz \
        ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10 \
        SLIDINGWINDOW:4:25 \
        MINLEN:45
```
   
`Trimmomatic` is a java program, so we invoke it by `java -jar $Trimmomatic`. In this case `$Trimmomatic` is an environmental variable set when you load the module that points to the actual jarfile for java to run. `PE` indicates we're running it in "paired-end" mode. `-threads 4` means we want to use 4 cpu threads to run it (make sure `-c 4` is set in the SLURM header!). Then we give the input and output fastq files. First the input R1 and R2, then their respective output files. The "singles" files are for reads whose mate pairs were too low quality to be retained. They are placed in a separate set of files. 

Finally, we specify the trimming parameters. `ILLUMINACLIP:/isg/shared/apps/Trimmomatic/0.36/adapters/TruSeq3-PE-2.fa:2:30:10` gives the adapter sequences we want to be trimmed off and a few parameters for the trimming algorithm. `SLIDINGWINDOW:4:25` means that if at any point in a scan of the sequence the average quality drops below 25 in a 4bp span, the rest of the read will be cut off. `MINLEN:45` indicates that the read should be dropped altogether if it is trimmed shorter than 45bp. For a more detailed explanation of these (and other) options, [see the documentation](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf.) 
At the end of the run, each run will produce **4** files, a *trimmed forward read file*, *trimmed reverse read file* and two *singles file*. Singles file will contain the reads whose pair failed to pass the thresholds set in the trimmommatic step. The following files will be produced at the end of the run:  

 
The summary of the reads will be in the `*.err` file, which will give how many reads is kept and how many have been discarded in each run.  
  
| Sample | Input records | Paired records kept | Forward surviving | Reverse surviving | Records dropped | Paired Kept (%) |   
| --- | --- | --- | --- | --- | --- | --- |   
| K21 | 19786941 | 15796097 | 2249710 | 577243 | 1163891 | 79.83% |
| K22 | 21277520 | 16946089 | 2365838 | 733446 | 1232147 | 79.64% |
| K23 | 22508598 | 18004494 | 2511606 | 751067 | 1241431 | 79.99% |
| K31 | 22058553 | 17905294 | 2212493 | 806168 | 1134598 | 81.17% |
| K32 | 20828110 | 16468780 | 2448484 | 759367 | 1151479 | 79.07% |
| K33 | 25096066 | 20321078 | 2601736 | 883807 | 1289445 | 80.97% |   
   


## 3. Assembling Transcriptomes
       
Our main goal here is to quantify gene expression in each of our samples. Because we used the illumina platform to do our sequencing, our read lengths (and fragment sizes) are much shorter than the average transcript. To identify which reads were derived from which genes, tally them up, and thus measure gene expression, the typical approach is to map the short reads back to a reference (either genome or transcriptome). Because we're working with a non-model organism, neither of those things are available. So here we'll assemble a transcriptome _de novo_. For a detailed review of genome assembly see [Simpson and Pop (2015)](https://www.annualreviews.org/doi/10.1146/annurev-genom-090314-050032). 

### De Novo Assembly using `Trinity`

To create a single reference transcriptome for all six samples, we'll first assemble each sample separately, pool all the resulting transcripts, cluster them into groups that hopefully represent single genes, and finally select a single "best" representative transcript for each gene. 

To do the assemblies we'll use the software [`Trinity`](https://github.com/trinityrnaseq/trinityrnaseq/wiki). Assembly requires a great deal of memory (RAM) and can take few days if the read set is large. The script to run `Trinity` is in the `03_Assembly` directory. You can run it from that directory by entering `sbatch trinity.sh` from the command line. Following is the command that we use to assemble each individual. Each of these assemblies will take between 6 and 15 hours. 
   
```bash
SAM=K21

Trinity --seqType fq \
  --left ../02_Quality_Control/trim_${SAM}_R1.fastq.gz \
  --right ../02_Quality_Control/trim_${SAM}_R2.fastq.gz \
  --min_contig_length 300 \
  --CPU 36 \
  --max_memory 100G \
  --output trinity_${SAM} \
  --full_cleanup  

```  

Trinity combines three independent software modules: Inchworm, Chrysalis, and Butterfly, applied sequentially to process large volumes of RNA-seq reads. Trinity partitions the sequence data into many individual de Bruijn graphs, each representing the transcriptional complexity at a given gene or locus, and then processes each graph independently to extract full-length splicing isoforms and to tease apart transcripts derived from paralogous genes. Per the [`Trinity` wiki](https://github.com/trinityrnaseq/trinityrnaseq/wiki), the process works like so:    
_Inchworm_ assembles the RNA-seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.   
   
_Chrysalis_ clusters the Inchworm contigs and constructs complete de Bruijn graphs for each cluster. Each cluster represents the full transcriptional complexity for a given gene (or sets of genes that share sequences in common). Chrysalis then partitions the full read set among these disjoint graphs.   
    
_Butterfly_ then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.   
   
During the `Trinity` run lots of files will be generated. These checkpoint files will help us to restart from that specific point if for some reason the program stops prematurely. Once the program ends successfully all these checkpoint files will be removed since we have requested a full cleanup using the `--full_cleanup` command. Clearing the files is very important as it will help us to remove all the unwanted files and also to keep the storage capacity and the number of files to a minimum. So at the end of a successful run we will end up with the following files: 
   
```
03_Assembly/
├── trinity_K21.Trinity.fasta
├── trinity_K21.Trinity.fasta.gene_trans_map
├── trinity_K22.Trinity.fasta
├── trinity_K22.Trinity.fasta.gene_trans_map
├── trinity_K23.Trinity.fasta
├── trinity_K23.Trinity.fasta.gene_trans_map
├── trinity_K31.Trinity.fasta
├── trinity_K31.Trinity.fasta.gene_trans_map
├── trinity_K32.Trinity.fasta
├── trinity_K32.Trinity.fasta.gene_trans_map
├── trinity_K33.Trinity.fasta
├── trinity_K33.Trinity.fasta.gene_trans_map
└── trinity.sh
```
   
So we will have six assemblies in fasta format, one for each sample.  
  
## 4. Identifying the Coding Regions   

Our goal in this portion of the tutorial is to create a single reference transcriptome containing one reference transcript per gene to quantify gene expression in all six samples against. We need to identify the coding regions of the transcripts for the first step in this process because we'll perform the clustering of redundant transcripts within and among samples using coding sequences. Later on, we'll also use amino acid sequences to do annotation using EnTAP. 
   
### Identifying coding regions using `TransDecoder` and `hmmer`  

Before looking for the coding regions we'll combine all the assemblies together. Because of `Trinity's` naming convention, there are likely to be identical sequence names appearing in each assembly. These redundant names would confuse `TransDecoder` and turn some of its output into nonsense. Worse, this error would occur silently and be propagated to the rest of our analyses. To deal with this, we'll simply add a sample ID prefix to each sequence name using the linux utility `sed`. For a single file:

```bash
SAM=K21
sed "s/>/>${SAM}_/g" ../03_Assembly/trinity_${SAM}.Trinity.fasta > ../03_Assembly/trinity_prefix_${SAM}.Trinity.fasta
```
Now we can concatenate the assemblies into a single file:

```bash
cat ../03_Assembly/trinity_prefix_* > ../03_Assembly/trinity_combine.fasta
``` 

Now that we have our reads assembled and combined together into the single file, we can run `TransDecoder` and `hmmer` to identify likely coding regions. There are three steps here. First, identify long [**open reading frames**]((https://en.wikipedia.org/wiki/Open_reading_frame)) (ORFs) using `TransDecoder.LongOrfs`. Many transcripts will have multiple ORFs that could correspond to a true coding sequence. Second, identify ORFs with homology to known proteins using [`hmmer`](http://hmmer.org/). Third, use `TransDecoder.Predict` to generate a final set of candidate coding regions. 


To find all long ORFs, we run `TransDecoder.LongOrfs` like this:
   
```
TransDecoder.LongOrfs -t ../Assembly/trinity_combine.fasta
```

By default it will identify ORFs that are at least 100 amino acids long (you can change this by using -m parameter). It will produce a folder called `trinity_combine.fasta.transdecoder_dir`.  

```
04_Coding_Regions
├── trinity_combine.fasta.transdecoder_dir
│   ├── base_freqs.dat
│   ├── longest_orfs.cds
│   ├── longest_orfs.gff3
│   └── longest_orfs.pep
```

Many transcripts are likely to have multiple ORFs, so it will be helpful to bring in evidence of homology to known proteins for the final predictions. We'll do this by searching the `Pfam` database. `Pfam` stands for "Protein families", and is a massive database with mountains of searchable information on, well, you guessed it, protein families. This will maximize the sensitivity for capturing the ORFs that have functional significance. We can scan the `Pfam` database using the software `hmmer`. The `Pfam` database is much too large to install on a local computer. We have a copy on Xanadu, however, in the directory `/isg/shared/databases/Pfam/Pfam-A.hmm`. We'll run `hmmer` like this:
   
```
hmmscan --cpu 16 \
        --domtblout pfam.domtblout \
        /isg/shared/databases/Pfam/Pfam-A.hmm \
        trinity_combine.fasta.transdecoder_dir/longest_orfs.pep
```

Once the run is completed it will create the following files in the directory.

```
04_Coding_Regions
├── pfam.domtblout
```
    
Lastly we use the 'TransDecoder.Predict' function along with our `hmmer` output to make final predictions about which ORFs in our transcripts are real.  
   
```bash
TransDecoder.Predict -t ../Assembly/trinity_combine.fasta \
        --retain_pfam_hits pfam.domtblout \
        --cpu 16
```   
   
This will add output to our `trinity_combine.fasta.transdecoder_dir`, which now looks like:
```
04_Coding_Regions/
├── pfam.domtblout
├── pipeliner.38689.cmds
├── pipeliner.5719.cmds
├── pipeliner.63894.cmds
├── trinity_combine.fasta.transdecoder.bed
├── trinity_combine.fasta.transdecoder.cds
├── trinity_combine.fasta.transdecoder_dir/
│   ├── base_freqs.dat
│   ├── hexamer.scores
│   ├── longest_orfs.cds
│   ├── longest_orfs.cds.best_candidates.gff3
│   ├── longest_orfs.cds.best_candidates.gff3.revised_starts.gff3
│   ├── longest_orfs.cds.scores
│   ├── longest_orfs.cds.top_500_longest
│   ├── longest_orfs.cds.top_longest_5000
│   ├── longest_orfs.cds.top_longest_5000.nr
│   ├── longest_orfs.gff3
│   └── longest_orfs.pep
├── trinity_combine.fasta.transdecoder.gff3
└── trinity_combine.fasta.transdecoder.pep
```

The full script is called [transdecoder.sh](/04_Coding_Regions/transdecoder.sh). It can be run from the `04_Coding_Regions` directory by entering `sbatch transdecoder.sh` on the command line. 
   

## 5. Determining and Removing Redundant Transcripts

De novo transcriptome assemblies are complex, containing both biological variation in the form of alternately spliced transcripts and nucleotide sequence variation, and technical issues such as fragmented transcripts. In this tutorial we have six de novo transcriptome assemblies. We're aiming to quantify gene expression at the gene level, so ideally we want to winnow out much of this complexity and create a single transcriptome with one transcript representing each underlying gene against which we can quantify gene expression for all six samples. In the previous step, we identified candidate coding regions for transcripts from all six samples. In this step, we'll cluster all those transcripts by amino acid sequence and select a single one to represent each cluster. 

### Clustering using `vsearch`

We will use `vsearch` to cluster the transcripts with similar sequences (similarity is set by the identity between the sequences --id) and then choose one representative transcript (the centroid). A more detailed account of the application can be found at: [vsearch](https://github.com/torognes/vsearch). The threshold for clustering in this example is set to 90% identity. We'll run `vsearch` like this: 
 
```bash
vsearch --threads 8 --log LOGFile \
        --cluster_fast ../Coding_Regions/trinity_combine.fasta.transdecoder.cds \
        --id 0.90 \
        --centroids centroids.fasta \
        --uc clusters.uc
```

The full script is called [vsearch.sh](/05_Clustering/vsearch.sh). It can be run from the `05_Clustering` directory by entering `sbatch vsearch.sh` on the command line. At the end of the run it will produce the following files:

```
05_Clustering/
├── centroids.fasta
├── clusters.uc
├── combine.fasta
└── LOGFile
```

The `centroids.fasta` file will contain the representative transcripts from the 6 assemblies. 
    

## 6. Evaluating the Assembly  

Now that we've settled on our reference transcriptome (`centroids.fasta`), we can assess its quality. We will benchmark it against a gene database using [rnaQUAST](http://cab.spbu.ru/software/rnaquast/). This will compare the transcripts with a gene database and will produce a summary report. 


```bash
rnaQUAST.py --transcripts ../05_Clustering/centroids.fasta \
	--gene_mark \
  --threads 8 \
  --output_dir Genemark
```

The full slurm script is called [rnaQuast.sh](/RNAQuast/rnaQuast.sh) and it can run from the `06_RNAQuast` directory by entering `sbatch rnaQuast.sh` on the command line. The resulting directory will look like this:

```
RNAQuast/
├── results/
│   ├── centroids_output/
│   │   ├── alignment_metrics.txt
│   │   └── txt_dir/
│   ├── logs/
│   │   └── rnaQUAST.log
│   ├── short_report.pdf
│   ├── short_report.tex
│   ├── short_report.tsv
│   └── short_report.txt
└── rnaQuast.sh
```

The program will produce various statistics and it will produce a short summary report as well as a long report, where you can look and assess the quality of your assembled transcriptome. Shown below are matrics which is produced in *alignment_metrics.txt* file. 

METRICS   |    centroids    
----  |   ---- 
Transcripts  |  64275   
Transcripts > 500 bp   |   28385  
Transcripts > 1000 bp   |  10593  
Average length of assembled transcripts  |   663.088  
Longest transcript   |   13230
Total length   |    42619973   
Transcript N50  |    407   


## 7. Functional Annotation

Now that we have our reference transcriptome, we want to annotate it. 

Here we'll use [`EnTAP`](https://entap.readthedocs.io/en/v0.9.0-beta/index.html). Our script to run `EnTAP` can be found in the directory `07_EnTAP` and can be run from that directory by entering `sbatch entap.sh` on the command line. 
  
The first step is to get the peptide sequences from our transdecoder output that correspond to the final reference transcriptome. We will use these as input for `EnTAP`. We can do this using `seqtk`:

```bash
grep -oP "(?<=>).*" $CENTROIDS >names.txt

CENTROIDS=../05_Clustering/centroids.fasta
PEPTIDES=../04_Coding_Regions/trinity_combine.fasta.transdecoder.pep

seqtk subseq $PEPTIDES names.txt >centroids.pep
``` 

The `subseq` module of `seqtk` extracts sequences from the first argument (a fasta formatted file) according to a list of names in the second argument. In this case we extract the the names from `centroids.fasta` using `grep` and a regular expression. 

The next step is to run `EnTAP`. To run `EnTAP` you need to have a configuration file set up in your working directory, where it points to the program paths. We have prepared a file (`entap_config.txt`) for you and it can be found in the `07_EnTAP/` directory.

The `EnTAP` command looks like this:

```bash
EnTAP --runP \
-i centroids.pep \
-d /isg/shared/databases/Diamond/RefSeq/plant.protein.faa.97.dmnd \
-d /isg/shared/databases/Diamond/Uniprot/uniprot_sprot.dmnd \
--ontology 0  \
--threads 8 \
-c bacteria \
-c fungi \
--taxon Larix
```

The `-d` flags specify the databases to use for the search. The `-c` flags will mark sequences that have hits to bacteria and fungi as possible contaminants, while the `--taxon` flag will favor annotations coming from congenerics in *Larix*. 
   
Once the job is done it will create a folder called “outfiles” which will contain the output of the program.   

```
EnTAP/
```

More information on EnTAP can be found in [EnTAP documentation](https://entap.readthedocs.io/en/v0.9.0-beta/index.html), which has a very comprehensive description.

   
## 8. Quantifying gene expression

The next step in our analysis is to quantify gene expression in each of our six samples. In this tutorial we'll use [`kallisto`](https://pachterlab.github.io/kallisto/), a very fast read pseudo-mapper. Without getting into detail, `kallisto` is fast, and referred to as a "pseudo-mapper" because it avoids the hard work of a base-level alignment of reads in favor of finding reads' approximate position in the reference transcriptome. This approximation is good enough for RNA-seq (but not for, e.g. variant calling). `kallisto` will produce estimated counts of reads mapping to transcripts, and estimated transcript abundances, normalized for transcript length. We'll use the estimated counts in later analyses. 

### Creating an index  

The first step is to index our reference transcriptome. An index allows possible read positions in the reference to be looked up very quickly. We'll be working in the `08_Counts` directory, and the index script `kallisto_index.sh` can be executed form that directory by entering `sbatch kallisto_index.sh` on the command line. 

```bash
kallisto index -i ../05_Clustering/centroids.fasta.index ../05_Clustering/centroids.fasta
```   

### Counting reads mapping to transcripts

Now we can use `kallisto` to quantify gene expression. Again from the `08_Counts` directory, we can run the script `kallisto_counts.sh` by entering `sbatch kallisto_counts.sh` on the command line. 

This script uses the `kallisto quant` module to run the quantification algorithm. Running it on a single sample looks like this:

```bash
SAM=K21

kallisto quant \
  -i ../05_Clustering/centroids.fasta.index \
  -o ${SAM} \
  -t 8 \
  ../02_Quality_Control/trim_${SAM}_R1.fastq.gz ../02_Quality_Control/trim_${SAM}_R2.fastq.gz
``` 
   
    
`kallisto` can process either paired-end or single-end reads. The default running mode is paired-end reads and requires a even number of FASTQ files, with pairs given as shown in above example. We provide the index with `-i` and the number of CPU threads with `-t`. 
  
The quantification algorithm will produce three output files: 
 *  abundance.h5  : HDF5 binary file   
	This contains run information, abundance estimates, bootstrap estimates and transcript lenght information
 *  abundance.tsv : plaintext file of the abundance estimates  
	This contains effective length, estimated counts and TPM values  
 *  run_info.json : information on the run  
  

Now if you look at the first few lines in the _abundance.tsv_ file using the `head` command. We can see that it has five columns: _geneID_, _gene length_, _effective gene length_, _estimated counts_ and _tpm_.  

```
target_id	length	eff_length	est_counts	tpm
TRINITY_DN27913_c0_g1_i1.p3	13230	13072.3	172	0.861807
TRINITY_DN27054_c1_g6_i1.p1	11508	11350.3	183.401	1.05834
TRINITY_DN26839_c0_g2_i1.p1	10935	10777.3	16.3293	0.099241
TRINITY_DN21012_c2_g1_i3.p1	10839	10681.3	172	1.05472
```

We can now take this information into our subsequent analyses. There we will look at our functional annotation and 

## 9. Differential expression analysis in R

The tutorial for this section is not yet written! But please see the extensive R script in `/09_R_analysis`. This script looks at transcript abundances across our taxonomic and functional annotations to identify the impact of contamination in this data set and finally does differential expression analysis with DESeq2 and gene ontology enrichment with goseq. 






 

 



