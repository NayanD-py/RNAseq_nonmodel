# THIS TUTORIAL IS UNDER CONSTRUCTION
  the scripts are broken after Step 8, quantifying gene expression

# RNA Sequence Analysis for Non Model Species Eastern larch (Tamarack)  
  
This repository is a usable, publicly available tutorial for analyzing differential expression data. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use `sbatch scriptname` after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use [this](https://bioinformatics.uconn.edu/unix-basics) handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](https://en.wikipedia.org/wiki/SAM_(file_format)), and [GFF3/GTF](https://en.wikipedia.org/wiki/General_feature_format). You can learn even more about each file format [here](https://bioinformatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/). If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one **[here](https://bioinformatics.uconn.edu/contact-us/)**.  
  
### Contents  
1. [Introduction](#1-introduction)  
2. [Quality Control](#2-quality-control)   
3. [Assembling Transcriptomes](#3-assembling-transcriptomes)  
4. [Identifying the Coding Regions](#4-identifying-the-coding-regions)  
5. [Determining and Removing Redundant Transcripts](#5-determining-and-removing-redundant-transcripts)
6. [Evaluating the Assembly](#6-evaluating-the-assembly)
7. [Functional Annotation](#7-functional-annotation) 
8. [Quantifying Gene Expression](#8-quantifying-gene-expression)
9. [Excluding contaminant genes](#9-excluding-contaminant-genes)
10. [Diffferential Expression](#10-diffferential-expression)  
       a.    [NOISeq](#a-differentially-expressed-genes-using-noiseq) 
       b.    [Gfold](#b-differentially-expressed-genes-using-gfold)   


## 1. Introduction  
  
The goal of this tutorial is to guide you through a differential gene expression analysis using RNA-seq data in a **non-model organism**. When an organism is called a **model** this usually means that very generous amounts of research have been performed on it, resulting in large pools of publicly available data. In bioinformatics this means there are pre-existing genomes or transcriptomes, gene annotations, and possibly a wealth of other useful information available in public databases. By contrast, when an organism is called **non-model** that often means those things are not available, and to do a genome-scale analysis, the researcher will have to produce them. 

So in this tutorial, in addition to the standard steps taken when those resources are available (e.g. in our [other RNA-seq tutorial](https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation)), we will also assemble and annotate a transcriptome, leading to a workflow that looks like this:

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

In this workflow we have separated each step into folders. In each folder is one or more script for each step along with the results after the scripts are run . When you clone the git repository, the following directory structure will be cloned into your working directory.   

So to follow the steps would be esay once you have cloned this git repository using the `clone` command:
```
git clone < git-repository.git >

```
Once you clone the repository you can see the following folder structure:  

```  
Eastern_larch/
├── 01_Raw_Reads
├── Quality_Control
├── Assembly
├── Coding_Regions
├── Clustering
├── RNAQuast
├── Index
├── Counts
├── Gfold 
├── NOISeq
└── EnTAP
```  
   
###  SLURM scripts   
This tutorial is tailored to the University of Connecticut's Xanadu computing cluster. Before beginning the tutorial, we'll need to understand a few aspects of the cluster. When first logging into Xanadu from your local terminal, you will be connected to a __submit node__. The submit node provides an interface for interacting with the cluster. Heavy computational work is not done on submit nodes, but rather on __compute nodes__. Never, under any circumstance do substantial computation on a submit node. You may inspect and move files around, or do light editing of scripts, but if you are running a command that uses more than a single processor, or takes more than a couple minutes you should be using a compute node. If you do use a submit node, your process could be killed, all of your work lost, and you may receive a nasty e-mail from a cluster administrator. 

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

The reads with which we will be working have been sequenced using [Illumina](https://www.illumina.com/techniques/sequencing.html). We'll assume here that you are familiar with the sequencing technology. Let's have a look at the content of one of our reads, which are in the `fastq` format. They are additionally compressed using `gzip`, so we'll use `zcat`, a pipe (`|`) and `head` to decompress and view the first sequence record:

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
## 8. Excluding Contaminant Genes



## 10. Diffferential Expression 

In this section we will demonstrate two methods of identifying differentially expressed genes: `Gfold` and `NOISeq`, but now that you have the count data, you could bring that to many other methods. To learn how to use `DESeq2`, for example, you could switch over to one of our [other tutorials](https://github.com/CBC-UCONN/RNA-seq-with-reference-genome-and-annotation#7-pairwise-differential-expression-with-counts-in-r-using-deseq2) and follow that from here.

Where as Gfold and NOISeq both can used to find differentially expressed genes when there are no-replicates as well as replicates avaliable. In this tutorial we will use Gfold, to find DE genes, when you do not have replicates and next we will use the NOISeq, show an example when you have replicates. 

### a. Differentially Expressed Genes using NOISeq  

Another program which is useful in finding differentially expressed genes when there are no replicates is the R package [NOISeq](https://bioconductor.org/packages/release/bioc/html/NOISeq.html). It can be used to get exploratory plots to evaluate count distribution, types of detected features, and differential expression between two conditions.   
  
In order to get the expressed genes using NOISeq program, you need to provide the count files as the input. For that we will use the kallisto generated count files which can be found in *abundance.tsv* files generated for each sample in the indexing step. As the abundance.tsv file will contain _GeneSymbol, GeneName, Read Count, Gene exon length and RPKM_ and we will use a *awk* command to grab only the _GeneSymbol_ and _Read Count_ from the *abundance.tsv* file.  

```
awk '{print $1 "\t" $4}' ../Counts/K23/abundance.tsv > K23.counts
awk '{print $1 "\t" $4}' ../Counts/K32/abundance.tsv > K32.counts
awk '{print $1 "\t" $4}' ../Counts/U32/abundance.tsv > U32.counts
awk '{print $1 "\t" $4}' ../Counts/U13/abundance.tsv > U13.counts
```   

We also need the transcript length information, for the NOISeq calculation. It can be gathered from the *abundance.tsv* file as well, using the following code:  

```awk
awk '{print $1 "\t" $2}' ../Counts/K23/abundance.tsv > length.txt
``` 

You can run this command in an interactive session or can run the [counts.sh](/NOISeq/counts.sh) script using shell in a interactive session.   
  
```bash
sh counts.sh
```   

This will produce counts files with two columns where it will contain the Gene Name and the counts associated with it. Also the length file with two columns which includes the transcript name and transcript length.  

```
NOISeq/
├── U13.counts
├── U32.counts
├── K23.counts
├── K32.counts
└── length.txt
```   

  
### NOISeq to produce differentially expressed genes   

This part can be run using the cluster or using your own PC or laptop. We recommend transferring the above counts files and the feature length file to your own computer and using the RStudio package.

To transfer the count files and the transcript feature length file from cluster to your location using a terminal window: 
```bash
scp <user name>@transfer.cam.uchc.edu:<PATH-to-NOISeq-Directory>/Eastern_larch/NOISeq/*.counts .
scp <user name>@transfer.cam.uchc.edu:<PATH-to-NOISeq-Directory>/Eastern_larch/NOISeq/length.txt .
```

**Prerequisites:** Before using the NOISeq package you need to make sure you have already downloaded the required packages including: NOISeq and dplyr. If not please download the necessary packages which are compatible with your version of R. 

When running the R code described below we will assume you have downloaded the count files to your R working directory.   

We will load the the *dplyr* library and will set the paths to the input count files and output files which will be generated during the execution of the program. The following script will write csv files and image files to the current directory, if you would like to direct it to a different location you are welcome to do so be changing the path of the output.  

```r
library(dplyr)

list.files()

#Set directory paths to working directory
count_dir <- getwd() # or appropriate path to the counts 
csv_out <- getwd() # or appropriate path to your output csv files
image_out <- getwd() # or appropriate path to your image out files
```   

Next we will create a dataframe to hold the count files by reading one at a time.  
```r
##  Read count files AND create a dataframe to hold the count data

#create a empty dataframe called m to merge the data into
m = data.frame()
# using for loop read all the count files in the count_dir path
for (i in list.files(pattern = ".counts")) {
  print(paste0("reading file: ", i))
  #read file as a data frame
  
  f <- read.table(i, sep = "\t", header = TRUE)
  #rename the columns
  colnames(f) <- c("gene_id", substr(i, 1, nchar(i)-7))
  #copy the data to another dataframe called f1
  f1 <- subset(f, select= c("gene_id", substr(i, 1, nchar(i)-7)))
  
  #if the m is empty just copy the f to m
  if(length(m) == 0){
    m = f1
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    m <- merge(m, f1, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f1)
}

#grab the rows from the 1st colum and use it as the row-names in the dataframe
rownames(m) <- m[,1]

# remove the column-1 (gene_ids) from the data frame using dplyr::select function
m <- select(m, "K23", "K32", "U13", "U32")
rm(f)
```   

Once the dataframe is created each column is represented by a sample and the rows with feature counts.  
 
Next we will prepare the metadata for the analysis. In here we can include the samples and different conditions we are hoping the evaluate. For this experiment we have two time points of data from a same location we will include that using TimePoint as a condition, which we would like to evaluate.   

```r
#################################################
## MetaData
#################################################
Sample = c("K32", "K23", "U13", "U32")
Condition = c("Killingworth", "Killingworth", "UConn", "UConn")
TimePoint = c("T2", "T3", "T3", "T2")
myfactors <- data.frame(Sample, Condition, TimePoint)
myfactors
```   

Once the myfactors dataframe is created it is important to check the columns in the dataframe(m) and mata-data rows are in the the same order.  

```r
# first check wheter all the columns in dataframe(m) and myfactor is present
all(colnames(m) %in% myfactors$Sample)

# Then check the order is correct
all(rownames(m) == myfactors$Sample[1])

# Then order them according to the Sample names
m <- m[, myfactors$Sample]

# Now check the order is correct after sorting 
all(colnames(m) == myfactors$Sample)
```   

We need to import the length information as a part of our meta data. 
```r
# import length information to a dataframe
df_length <- read.table("length.txt", sep = "\t", header = TRUE, row.names = 1)

# create a vector to hold the length information
mylength <- setNames(object = df_length$length, row.names(df_length))
head(mylength)
```
  

As the next step we will create the NOISeq object using the count dataframe and factors dataframe created above.  

```r
#####################################################
##  NOISeq          
#####################################################

library(NOISeq)

# Creating a NOISeq object
mydata <- readData(data = m, factors = myfactors, length = mylength)
mydata

head(assayData(mydata)$exprs)
head(featureData(mydata)@data)
head(pData(mydata))
```   

#### Quality Control of Count Data  
Before proceeding with the analysis of data it is important to see the quality of the experimental data. This includes determining possible biases of the data. In this tutorial we have shown the following plots to identify the quality of our count data.   
1.   Count distribution 
2.   Lengthbiases   
3.   Batch effects   

Before generating any type of plots, *dat* function must be applied on the input data which is the NOISeq object to obtain the data which is needed. This can be passed using the _type_ argument. Once the data is generated to plot, image can be generated using the _explo.plot_ function.

##### 1. Count distribution plot  
As a part of a quality control step before procedding with the analysis it is important to visualize the possible bias of the sample. As a part of this we can detect the count distribution of samples.   

```r
#####################################################
## Count Distribution 
#####################################################
countsbio = dat(mydata, factor = NULL, type = "countsbio")
explo.plot(countsbio, toplot = 1, samples = NULL, plottype = "boxplot")
```
![](images/Count_distribution_plot.png)

##### 2. Length bias plot  
Lengthbias plot will discribe the relation of the feature length towards the expression values. Once caluclated, both model p-value and coefficient of determination (R2) are shown in the plot with the regression curve. 

```r
##############
## Length Bias
##############
mylenghtbias = dat(mydata, factor = "TimePoint", type = "lengthbias")
explo.plot(mylenghtbias, samples = NULL, toplot = "global")
```

![](images/Lengthbias_plot.png)   

Significant p-value with a high R2 value will indicate a expression is dependent on the feature length and the curve shows its dependence.  


##### 3.  PCA plots 
Principle component analysis (PCA) plots can be used to visualize how the samples are clustered acoording to the experimental design.   


```r
##############
## PCA
##############
myPCA = dat(mydata, type = "PCA")
explo.plot(myPCA, factor= "TimePoint")
```   

![](images/PCA_plot_Timepoint.png)  


#### NOISeq-bio  

NOISeq method can be used to compute differential expression on data set with technical replicates (**NOISeq-real**) or without replicates (**NOISeq-sim**). Also NOISeq method can be applied when there are biological replicates (**NOISeq-bio**).   
  
In here we will show you how to calculate the differential expression, when there is biological replicates. We will take into account the TimePoint as our factor where it has biological replicates at T2 and T3.   

```r  
###########################
## Differential Expression
###########################
# lc = 1 RPKM values with length correction
# lc = 0 no length correction applied
mynoiseq.bio <- noiseqbio(mydata,
                          factor = "TimePoint",
                          norm = "rpkm",
                          random.seed = 12345,
                          lc = 1)

head(mynoiseq.bio@results[[1]])
```  

**noiseqbio** function can be used to calculate differential expression, which have biological replicates.  `norm` is the normalization method to be used, which can be rpkm by default, or uqua (upper quartile), tmm (trimmed mean of M) or n which is no normalization. In here `lc = 1 ` is used to calculate the RPKM values with the length correction.    

This will produce a noiseq object which contains following elements:

```bash
                              T2_mean    T3_mean      theta      prob     log2FC length
TRINITY_DN0_c0_g1_i1.p1     14.014433 10.6701081  0.0880846 0.5045883  0.3933386    883
TRINITY_DN1_c0_g1_i1.p2      9.889114 29.7867351 -0.3592379 0.9703946 -1.5907568    711
TRINITY_DN1000_c0_g1_i1.p1   1.983390  0.9570247  0.1529210 0.6428457  1.0513405    714
TRINITY_DN10000_c0_g1_i2.p1 14.294755 19.6743936 -0.1065910 0.7681836 -0.4608333   1050
TRINITY_DN10001_c0_g1_i1.p2 84.469598  3.6968466  0.9136529 0.9971911  4.5140651    458
TRINITY_DN10001_c0_g1_i2.p1 43.963212 53.6521530 -0.1153515 0.7883396 -0.2873393   1398
```   

#### Selecting differentially expressed features  

The next step would be, how to select the differentially expressed features. This can be done using **degenes** function and applying a threshold using the q value. With the argument M, we can choose if we want all the differentially expressed genes (NULL), only the differentially expressed features that are more expressed in condition 1 than in condition 2 (M="up") or only the features which are under expressed in condition 1 with regard to condition 2 (M = "down").  

```r
mynoiseq.bio.deg = degenes(mynoiseq.bio, q = 0.9 , M= NULL)
```  
```
[1] "47586 differentially expressed features"
```

```r 
mynoiseq.bio.deg_up = degenes(mynoiseq.bio, q = 0.9 , M= "up")
head(mynoiseq.bio.deg_up)
```  
```
[1] "23355 differentially expressed features (up in first condition)"
                               T2_mean    T3_mean    theta prob    log2FC length
TRINITY_DN57434_c1_g3_i1.p1 2899.89293 0.08799537 2.785696    1 15.008213    447
TRINITY_DN58863_c0_g2_i9.p1   24.55380 0.02478509 1.712863    1  9.952258   1587
TRINITY_DN57127_c0_g1_i1.p1   52.09588 0.08043748 1.713380    1  9.339086    489
TRINITY_DN50342_c0_g1_i1.p1   53.81552 0.08404686 1.713626    1  9.322613    468
TRINITY_DN56349_c0_g1_i1.p1   42.69485 0.05853263 1.713635    1  9.510605    672
TRINITY_DN56871_c1_g1_i3.p2   52.80055 0.08143671 1.714048    1  9.340658    483
```  


```r
mynoiseq.bio.deg_down = degenes(mynoiseq.bio, q = 0.9 , M= "down")
head(mynoiseq.bio.deg_down)
```  
```
[1] "24231 differentially expressed features (down in first condition)"
                                  T2_mean     T3_mean     theta      prob      log2FC length
TRINITY_DN27469_c0_g1_i1.p1  6.879962e+02 1359.342646 -3.810915 1.0000000  -0.9824367    412
TRINITY_DN13350_c0_g2_i1.p1  1.207381e-07    4.437976 -3.351256 1.0000000 -25.1315171    327
TRINITY_DN21734_c1_g1_i5.p1  1.811139e-06   15.734521 -3.223580 1.0000000 -23.0505324    663
TRINITY_DN20971_c0_g1_i26.p1 3.058669e+02  643.731754 -2.358441 1.0000000  -1.0735556    372
TRINITY_DN25591_c0_g1_i1.p1  1.680417e-01  132.087619 -1.809176 0.9997979  -9.6184604    384
TRINITY_DN27792_c0_g1_i23.p1 6.759308e-01   40.725362 -1.798873 0.9987484  -5.9129081    675
```  

These features can be written into a csv file using the following command.
```r
prefix = "T2_T3_noiseq"
write.csv(mynoiseq.bio.deg_up, file = paste0(csv_out, "/" ,prefix, "_DEgenes_up.csv"))
write.csv(mynoiseq.bio.deg_down, file = paste0(csv_out, "/", prefix, "_DEgenes_down.csv"))
``` 

### b. Differentially Expressed Genes using Gfold   

   
Here we are trying to get differentially expressed genes between two conditions with only a single replicate for each condition. In such situations [Gfold](https://zhanglab.tongji.edu.cn/softwares/GFOLD/index.html) is very useful. Gfold uses a Bayesian model to rank genes according to estimated fold change in expression in across treatments. The Gfold estimate of fold change balances the size of the fold change against its significance, reducing the noise from genes with low read counts. 
 
In order to get fold change using Gfold program, you need to provide a count file in a particular format. For the Gfold program the count file should contain 5 columns, where it should contain _GeneSymbol_, _GeneName_, _Read Count_, _Gene exon length_ and _RPKM_. In this input most important columns are _Gene Symbol_ and _Read Count_ information.   
  
We will use the kallisto generated _abundance.tsv_ file and reformat it, so it will contain the fields which Gfold program needs as its input. For reformating you can use any programing language that you like. But in here we will be using `awk` to manipulate these columns. AWK is a very powerful language which allows you to get useful information from text files.  
   
```awk
awk '{print $1 "\t" $1 "\t" $4 "\t" $2 "\t" $5 }' ../Counts/K23/abundance.tsv > K23.read_cnt
awk '{print $1 "\t" $1 "\t" $4 "\t" $2 "\t" $5 }' ../Counts/K32/abundance.tsv > K32.read_cnt
```  

You can run this command in a interative session or can run the [count2gfoldCounts.sh](/Gfold/count2gfoldCounts.sh) script using shell in a interative session. 
```bash
sh count2gfoldCounts.sh
```  

This will produce counts files which are in Gfold format.  
```
Gfold/
├── K23.read_cnt
└── K32.read_cnt
```   
   
Gfold program does not take the header row so either you have to delete the first row or comment out the header before you run the program. 
Now if you look at any of these files it will now contain five columns as follows:
```
TRINITY_DN27913_c0_g1_i1.p3	TRINITY_DN27913_c0_g1_i1.p3	3561	13230	12.3338
TRINITY_DN27054_c1_g6_i1.p1	TRINITY_DN27054_c1_g6_i1.p1	3895.81	11508	15.5417
TRINITY_DN26839_c0_g2_i1.p1	TRINITY_DN26839_c0_g2_i1.p1	1220.95	10935	5.12987
TRINITY_DN21012_c2_g1_i3.p1	TRINITY_DN21012_c2_g1_i3.p1	3349	10839	14.1975
TRINITY_DN17708_c0_g1_i3.p1	TRINITY_DN17708_c0_g1_i3.p1	967	9297	4.79158
```

### Using Gfold get differentially expressed genes  

Now since we have the counts files in the correct format, we will run the Gfold using the following command. We will be using the K23.read_cnt and K32.read_cnt as input files.
```bash
module load gfold/1.1.4

gfold diff -s1 K32 -s2 K23 -suf .read_cnt -o K32_vs_K23.diff
``` 

Usage information of the gfold:
```
gfold diff [options]

-s1      sample-1
-s2      sample-2
-o	 output file name
-suf	 input files extention

OPTIONAL FLAGS:
-acc <T/F>	When no replicate is available, whether to use accurate method to calculate GFOLD value. T stands for 	accurate which depends on sequencing depth and slower, F stands for MCMC. Default T. For job diff only.

-sc <num>	The significant cutoff for fold change. Default 0.01. For job diff only.

-norm <Count/DESeq/NO>  The way to do normalization. 'Count' stands for normalization by total number of mapped reads. 'DESeq' stands for the normalization proposed by DESeq. 'NO' stands for no normalization. You can also specifiy a list of normalization constant separated by commas. E.g. 1.2,2.1,1.0,2.0. Note that the number of constants should be the same as the total number of samples (group1 and group2) and the order should be for -s1 followed by for -s2. GFOLD using normalization constants not by directly multiplication (scaling up) nor division (scaling down). The normalization constants will be built into the model. In the model, division or multiplication has no difference. Default 'DESeq'.
```
   
The complete slurm script is called [gfold.sh](/Gfold/gfold.sh) which is stored in the **Gfold/** directory. Running Gfold will generate the following files, which contains the fold change value between the two conditions.
```
Gfold/
├── K32_vs_K23.diff
└── K32_vs_K23.diff.ext
```   
   
The Gfold value can be considered a log2-fold change value, where positive/negative value will indicate whether a gene is up/down regulated. The first few lines in the *K32_vs_K23.diff* will be like:  

```
#GeneSymbol     GeneName        GFOLD(0.01)     E-FDR   log2fdc 1stRPKM 2ndRPKM
TRINITY_DN27913_c0_g1_i1.p3     TRINITY_DN27913_c0_g1_i1.p3     2.49054 1       2.74244 2.33862 22.549
TRINITY_DN27054_c1_g6_i1.p1     TRINITY_DN27054_c1_g6_i1.p1     2.53847 1       2.78281 2.8605  28.3545
TRINITY_DN26839_c0_g2_i1.p1     TRINITY_DN26839_c0_g2_i1.p1     3.81393 1       4.54498 0.263204        9.34665
TRINITY_DN21012_c2_g1_i3.p1     TRINITY_DN21012_c2_g1_i3.p1     2.4016  1       2.65392 2.8545  25.8846
TRINITY_DN17708_c0_g1_i3.p1     TRINITY_DN17708_c0_g1_i3.p1     2.27485 1       2.74287 0.890033        8.71362
```  
The GFOLD output has 7 coloumns 
1. **GeneSymbol** Information on gene symbol.
2. **GeneName** Information on gene name.
3. **GFOLD** column provides fold change values in log2 format and can be used to obtain a biological meaningful ranking of genes. Any gene that passes the significance cutoff (p-value) of 0.01 and shows 2 or more fold change in expression have values indicated against them.  Genes not satisfying these 2 criterias have value=0 against them. Values are calculated as log2(s2/s1).  The significance cut off can be set by using -sc flag. Since the values are in log2 format, the cut off starts at +1 for upregulated genes and -1 for downregulated genes.  Genes with values greater than +1 (1.2435, 2.4982, 3.53474 etc) have 2 fold increase in expression in s2 samples wheres values less than -1 (-1.6584, -2.0078, -4.6768 etc) will have 2 fold lower or lesser expression in s2 compared to s1.
4. **E-FDR** column represents the FDR values calculated to correct for multiple testing.  In absence of replicates the value is set to 1 as seen above.  This column will have other values if we have replicates in our study.
5. **log2fc** column have log2 fold change obtained from s2/s1 for all genes even for those who doesnot pass the significance cut off and have lower or no change in expression between conditions. These values are slightly different from GFOLD because fold change is based on the sampled expression level from the posterior distribution bt taking in account gene lengths.
6. **1-RPKM** represent RPKM values for genes in s1 sample.
7. **2-RPKM** represent RPKM values for genes in s2 sample.







 

 



