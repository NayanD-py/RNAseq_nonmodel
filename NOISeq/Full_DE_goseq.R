#################################################
# Filtering, DE, and GO enrichment analysis
#################################################

# The purpose of this script is to read in the count data and
# the transcript annotations, filter out transcripts that are
# likely to be contaminants, and do a DE and GO enrichment. 
# We'll use NOISeq and DESeq2, and write out data to be analyzed
# using GFOLD. 


#################################################
# Load necessary libraries
#################################################

library(tidyverse)
library(NOISeq)
library(DESeq2)
library(pheatmap)
library(ggrepel)


#################################################
# Set some directories for input/output files
#################################################

list.files()

#Set directory paths to working directory
count_dir <- "08_Counts" # or appropriate path to the counts 


#################################################
# Read in count data
#################################################

# create a empty dataframe called co to merge the data into
co <- data.frame()

# using for loop read all the count files in the count_dir path
# iterate over each "abundance.tsv" file
# merge each file, i, into a single data frame

count_files <- list.files(path = count_dir,recursive=TRUE,pattern="abundance.tsv",full.names=TRUE)

for (i in count_files) {

  # print the file that is being loaded
  print(paste0("reading file: ", i))

  # read file i as a data frame
  f <- read.table(i, sep = "\t", header = TRUE)[,c(1,4)]

  # rename the columns
  colnames(f) <- c("gene_id", substr(i, 11, 13))

  #if the counts is empty just copy the f to m
  if(length(co) == 0){
    co <- f
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    co <- merge(co, f, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f)
}

#grab the rows from the 1st column and use it as the row-names in the dataframe
rownames(co) <- co[,1]

# create a vector containing transcript lengths
# read in length
df_length <- read.table(count_files[1], sep = "\t", header = TRUE)[,1:2]

# create named vector
mylength <- df_length$length
names(mylength) <- df_length[,1]


#################################################
# Read in annotation data
#################################################

# Here we read in the table containing our gene annotations from EnTAP
  # the sep, quote, and comment.char options ensure this complex table is parsed correctly
ann <- read.table("final_annotations_lvl3.tsv",header=TRUE, sep="\t", quote="", comment.char="")

# merge the counts data frame with the annotations
tab <- merge(co, ann, by.x = "gene_id", by.y = "Query.Sequence")

# EnTAP flags possible contaminant transcripts, but leaves the field blank if
# it does not find a hit for the transcript, change the blank to "U" for "Unknown"
tab[tab$Contaminant=="","Contaminant"] <- "U"

# Do the same for taxonomic scope
tab[tab$Tax.Scope=="","Tax.Scope"] <- "Unknown"

# add a column, Prevalence for the number of samples a gene is expressed in
tab$Prevalence <- rowSums(tab[,2:7]>0)


#################################################
# Summarize expression levels by contaminant status and taxonomy
#################################################

# get a table summarizing how many RNA fragments map by contaminant status
bycontam_wide <- group_by(tab,Contaminant) %>% summarize(across(K21:K33,sum)) %>% data.frame()
# pivot the table so that it can be easily turned into a barplot
bycontam <- bycontam_wide %>% pivot_longer(!Contaminant,names_to="Sample",values_to="Count")
# order the factors
bycontam$Contaminant <- factor(bycontam$Contaminant,levels=c("U","Yes","No"))

# stacked barplot of contaminant status for each sample, absolute values
ggplot(bycontam,aes(fill=Contaminant,y=Count,x=Sample)) +
  geom_bar(position="stack",stat="identity")

# stacked barplot of contaminant status for each sample, proprortions
ggplot(bycontam,aes(fill=Contaminant,y=Count,x=Sample)) +
  geom_bar(position="fill",stat="identity")

# grouped barplot of prevalence frequency by contaminant status
  # must call "dplyr::count" because tidyverse count was masked by another package

# extract prevalence frequencies by contaminant status
prevalence <- data.frame(Contaminant=tab$Contaminant,Prevalence=tab$Prevalence) %>% 
  group_by(Contaminant) %>% 
  dplyr::count(Prevalence) %>% 
  data.frame()

# grouped barplot
ggplot(prevalence,aes(fill=Contaminant,y=n,x=Prevalence)) +
  geom_bar(position="dodge",stat="identity")


# get a table summarizing how many RNA fragments map by taxonomic scope
bytax_wide <- group_by(tab,Tax.Scope) %>% summarize(across(K21:K33,sum)) %>% data.frame()
# pivot the table so that it can be easily turned into a barplot
bytax <- bytax_wide %>% pivot_longer(!Tax.Scope,names_to="Sample",values_to="Count")

# stacked barplot of taxonomic scope status for each sample, absolute values
ggplot(bytax,aes(fill=Tax.Scope,y=Count,x=Sample)) +
  geom_bar(position="stack",stat="identity")

# extract prevalence frequencies by taxonomic scope
prevalence <- data.frame(Tax.Scope=tab$Tax.Scope,Prevalence=tab$Prevalence) %>% 
  group_by(Tax.Scope) %>% 
  dplyr::count(Prevalence) %>% 
  data.frame()

# grouped barplot
ggplot(prevalence,aes(fill=Tax.Scope,y=n,x=Prevalence)) +
  geom_bar(position="dodge",stat="identity")


#################################################
# Use only transcripts flagged as Viridiplantae (exclude Unknowns)
#################################################

# create data frame 'm' where columns are counts, rownames are transcript IDs
m <- filter(tab,Tax.Scope=="Viridiplantae")[,1:7]
rownames(m) <- m[,1]
m <- m[,-1]

# also subset transcript lengths
mylength <- mylength[rownames(m)]


#################################################
# Create metadata table
#################################################

# here we create the table
Sample = c("K21","K22","K23","K31","K32","K33")
TimePoint = c("T2","T2","T2","T3","T3","T3")
myfactors <- data.frame(Sample, TimePoint)
myfactors

# the order and identity of samples in data frame m and myfactor must be the same
# ensure all the columns in dataframe m and myfactor is present
all(colnames(m) %in% myfactors$Sample)

# Then check the order is the same in both objects
all(colnames(m) == myfactors$Sample)

# If are not, so order them according to the sample names
m <- m[, myfactors$Sample]

# Now check the order is correct after sorting 
all(colnames(m) == myfactors$Sample)


#################################################
# Preliminary data exploration and filtering
#################################################

# how many transcripts did we quantify expression for?
dim(m)

# how many (non-contaminant) transcripts have expression for n samples?
prevalence <- filter(tab,Tax.Scope=="Viridiplantae")[,"Prevalence"]
table(prevalence)

# what proportion of expression data maps to low prevalence transcripts?
colSums(m[prevalence<2,])/colSums(m)

# what is the distribution of total expression across samples?
rowSums(m) %>% log(.,10) %>% hist(.,100)

# Before moving forward, we're going to exclude genes that are most likely to be useless in the analysis
# let's throw out all the transcripts with
  # a) less than 20 fragments mapping in total across all samples
  # b) that have non-zero expression for only a single sample

subg <- prevalence > 1 & rowSums(m) > 20
table(subg)

m <- m[subg,]
mylength <- mylength[subg]


#################################################
# First we'll analyze our data using DESeq2
#################################################

# we've already created the necessary component objects
  # use them to create a "DESeqDataSet" object
  # we are rounding 'm' because DESeq2 expects integer counts,
  # but Kallisto estimates the counts, resulting in fractional numbers

ddsHTSeq <- DESeqDataSetFromMatrix(
  countData = round(m),
  colData = myfactors,
  design = ~ TimePoint
  )

######################################################
# Reset treatment factors
######################################################

# To see the levels as they are now:
ddsHTSeq$TimePoint

# To replace the order with one of your choosing, create a vector with the order you want:
treatments <- c("T2","T3")

# Then reset the factor levels:
ddsHTSeq$TimePoint <- factor(ddsHTSeq$TimePoint, levels = treatments)

# verify the order
ddsHTSeq$TimePoint

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsHTSeq)

######################################################
# Get a table of results
######################################################

# get results table
resDE <- results(dds)

# get a quick summary of the table
summary(resDE)

# check out the first few lines
head(resDE)



######################################################
# Get a table of shrunken log2 fold changes
######################################################

# see coefficient names:
resultsNames(dds)

# get shrunken log fold changes, specifying the coefficient 
res_shrink <- lfcShrink(dds, coef="TimePoint_T3_vs_T2", type="ashr")

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=resDE$log2FoldChange,
  y=res_shrink$log2FoldChange,pch=20,
  cex=.2,
  col=1+(resDE$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change",
  xlim=c(-10,15),
  ylim=c(-10,15)
  )
abline(0,1)

points(
  x=resDE$log2FoldChange,
  y=res_shrink$log2FoldChange,pch=20,
  cex=.2,
  col=4 * is.na(resDE$padj)
  )

# get the top 20 genes by shrunken log2 fold change
top20 <- order(-abs(res_shrink$log2FoldChange))[1:20]
res_shrink[top20,]


######################################################
# Data visualization
######################################################

# MA plot
plotMA(res_shrink, ylim=c(-4,4))

# distribution of log2 fold changes:
  # there should be a peak at 0
hist(resDE$log2FoldChange,breaks=200)
abline(v=0,col="red",lwd=2)

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.2,
     col=(log_padj > 3)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="TimePoint")

# alternatively, using ggplot

dat <- plotPCA(vsd, intgroup="TimePoint",returnData=TRUE)

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group))
p <- p + geom_point() + geom_text_repel(aes(label=name))
p

##############

# heatmap of DE genes

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# get top 50 log fold change genes (excluding cook's cutoff outliers)
top50 <- data.frame(res_shrink) %>%
  filter(!is.na(padj)) %>% 
  arrange(-log2FoldChange) %>% 
  rownames() %>% 
  head(.,n=50)

df <- data.frame(colData(dds)[,"TimePoint"])
  rownames(df) <- colnames(dds)
  colnames(df) <- "TimePoint"

pheatmap(
  assay(rld)[top50,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
  )

##############

# plot counts for individual genes
  # get an ordered list of genes
  # exclude cook's cutoff genes
l2fc_ord <- data.frame(res_shrink) %>%
  filter(!is.na(padj)) %>% 
  arrange(-log2FoldChange) %>% 
  rownames()

plotCounts(dds, gene=l2fc_ord[1], intgroup="TimePoint")





#################################################
# Next we'll analyze our data using NOISeq
#################################################




#################################################
# Create the NOISeq data object
#################################################

# Creating a NOISeq object
mydata <- readData(data = round(m), factors = myfactors, length = mylength)
mydata

# access expression values
head(assayData(mydata)$exprs)
# access length values
head(featureData(mydata)@data)
# access sample ID to treatment mapping
head(pData(mydata))


#####################################################
## Exploratory data analysis
#####################################################

# many of these are similar to ones we can create using DESeq2

# plot count distributions by sample:
countsbio <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(countsbio, toplot = 1, samples = NULL, plottype = "boxplot")

# do we need to apply among-sample normalization (the answer is almost always yes!)
mycd <- dat(mydata, type = "cd", norm = FALSE, refColumn = 1)

# is there length bias in expression values?
  # note that points are averages across many genes
  # length does not really explain 80-90% of count variation
mylengthbias <- dat(mydata, factor = "TimePoint", type = "lengthbias")
explo.plot(mylengthbias, samples = NULL, toplot = "global")

# plot a PCA
  # what do you make of PC1?!
myPCA <- dat(mydata, type = "PCA", logtransf=FALSE)
explo.plot(myPCA, factor= "TimePoint")

#####################################################
## Differential Expression
#####################################################
# lc = 1 RPKM corrects for fragment length
# lc = 0 no length correction applied
# norm = "tmm" uses TMM normalization to normalize across samples

mynoiseq.bio <- noiseqbio(mydata, factor = "TimePoint", norm = "tmm", lc = 1)

resNOI <- mynoiseq.bio@results[[1]]
head(resNOI)

# MA plot

logmean <- log(rowSums(resNOI[,1:2]),10)
logFC <- resNOI$log2FC
sig <- resNOI$prob > 0.999
plot(logmean,logFC,pch=20,cex=.3,col=sig+1)

hist(resNOI[,]$log2FC,breaks=200)
abline(v=0,lwd=2,col="red")

# let's make a heatmap of some DE genes:

# get log-transformed (+1) normalized counts
myTMM <- (tmm(assayData(mydata)$exprs, long = 1000, lc = 1)+1) %>% log(.,10)

# get 100 genes with prob(DE)=1
top100 <- which(sig)[1:200]
df <- data.frame(myfactors[,-1])
rownames(df) <- colnames(myTMM)
#colnames(df) <- "condition"
pheatmap(
  myTMM[top100,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df,
  fontsize=5
)

# multidimensional scaling plot (similar to PCA)
mds <- t(myTMM) %>% dist() %>% cmdscale()
mds <- cbind(myfactors,mds)
colnames(mds)[3:4] <- c("mds1","mds2")

p <- ggplot(mds, aes(mds1, mds2, label = Sample)) +
  geom_point(aes(color = factor(TimePoint)))

p2 <- p + geom_text_repel(size=2) + labs(title = "mds plot of expression data")
p2


plot(log(res$T2_mean+1,2),log(res$T3_mean+1,2),pch=20,cex=.2,col=rgb(0,0,0,.4))
abline(0,1)

plot(myTMM[,3],myTMM[,4],pch=20,cex=.2,col=rgb(0,0,0,.4))


#######
# let's compare the results
#######

dat <- data.frame(NOISeq=-log(1-resNOI$prob,10),DESeq2=-log(resDE$padj,10),color=(as.numeric(resNOI$prob>0.99) + 2*(resDE$padj < 0.1))+1)
dat[,3] <- c("ns","NOI_sig","DE_sig","NOI & DE")[dat[,3]]

table(dat[,3])

ggplot(dat,aes(color=color,x=NOISeq,y=DESeq2)) +
  geom_point()



###############################################
# gene ontology enrichment analysis 
###############################################

library(goseq)

# get annotations (this is only a subset for this demonstration)
ann <- read.table("final_annotations_no_contam_lvl0.tsv",sep="\t",quote="",header=TRUE)

# get gfold results
gfold <- read.table("U32_vs_U13.diff",skip=1)
colnames(gfold) <- c("Query.Sequence","Gene.Name","gfold","eFDR","log2FC","RPKM_1","RPKM_2")

# combine these into a single data frame using only transcripts in common to both
tab <- inner_join(gfold,ann,by="Query.Sequence")

# select only transcripts for which an annotation exists (here we'll use a non-missing blast e-value as a shorthand)
# this will be our total gene set within which we'll look for GO terms associated with genes with extreme gfold scores
# note that here we only have 6810 transcripts, but for an actual analysis you will most likely have more
tab <- tab[!is.na(tab$E.Value),]

# now we need to format the data for goseq
  # a vector of 1/0 for each gene, indicating DE/not DE
  # a vector of transcript lengths (the method tries to account for this source of bias)
  # a table of transcript ID to category IDs (in this case GO term IDs) 

# DE/not DE vector: in this case let's select genes with gfold scores > 4 or < -4
  # 576 genes
de <- as.numeric(tab$gfold > 4 | tab$gfold < -4)
names(de) <- tab[,1]

# vector of transcript lengths
len <- mylength[tab[,1]]

# mapping of transcript IDs to category IDs

# a function to parse a single set of entap go terms
# a new version of entap will output the annotation pre-formatted like this
parse_entap_goterms <- function(goterms,geneid){
  if(goterms==""){return(c(geneid,NA,NA))}
  str_split(goterms,regex(",(?=GO)"))[[1]] %>%
    str_split(.,regex("(?<=GO:[0-9]{7})-")) %>% 
    do.call(rbind,.) %>% 
    cbind(geneid,.)
  }

# loop through each transcript and collect the GO terms
goterms <- c()
for(i in 1:dim(tab)[1]){
  
  x <- parse_entap_goterms(tab$GO.Biological[i], tab[i,1])
  y <- parse_entap_goterms(tab$GO.Molecular[i], tab[i,1])
  z <- parse_entap_goterms(tab$GO.Cellular[i], tab[i,1])
  
  goterms <- rbind(goterms,x,y,z)
  if((i %% 100)==0){print(i)}
  }

# extract and bind levels
goterms <- data.frame(goterms,levels=str_extract(goterms[,3],regex("(?<=\\(L=)[0-9]")))

# select level 3 GO terms, get only the transcript ID and go term ID
go <- data.frame(goterms[goterms$levels=="3",1:2])

# now we can start the analysis
# first try to account for transcript length bias by calculating the
# probability of being DE based purely on gene length
pwf <- nullp(DEgenes=de,bias.data=len)

GO.wall <- goseq(pwf=pwf,gene2cat=go)

# do FDR correction on p-values using Benjamini-Hochberg, add to output object
GO.wall <- cbind(
  GO.wall,
  padj_overrepresented=p.adjust(GO.wall$over_represented_pvalue, method="BH"),
  padj_underrepresented=p.adjust(GO.wall$under_represented_pvalue, method="BH")
)

# explore the results

head(GO.wall)


# identify transcript IDs associated with the top enriched GO term
g <- which(go$V2==GO.wall[1,1])
gids <- go[g,1]

# get their gfold results
gfold[gfold[,1] %in% gids,]

# get gene descriptions
tab[tab[,1] %in% gids, c("Query.Sequence","Description")]

# plot gfold scores for those genes, sorted
ord <- order(gfold[gfold[,1] %in% gids,]$gfold)
plot(gfold[gfold[,1] %in% gids,]$gfold[ord],
     ylab="gfold scores of genes in top enriched GO term",
     pch=20,cex=.5)
abline(h=0,lwd=2,lty=2,col="gray")




tab2 <- tab %>% 
  mutate(GO.Biological = str_split(GO.Biological,regex(".(?=GO:)"))) %>%   
  unnest(GO.Biological)

