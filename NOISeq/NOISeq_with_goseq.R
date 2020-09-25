#################################################
# Load necessary libraries
#################################################

library(tidyverse)
library(NOISeq)
library(pheatmap)
library(ggrepel)

#################################################
# Set some directories for input/output files
#################################################

list.files()

#Set directory paths to working directory
count_dir <- getwd() # or appropriate path to the counts 
csv_out <- getwd() # or appropriate path to your output csv files
image_out <- getwd() # or appropriate path to your image out files

#################################################
# Read in count data
#################################################

# create a empty dataframe called m to merge the data into
m <- data.frame()

# using for loop read all the count files in the count_dir path
# iterate over each ".counts" file
# merge each file, i, into a single data frame
for (i in list.files(path = count_dir, pattern = ".counts")) {

  # print the file that is being loaded
  print(paste0("reading file: ", i))

  # read file i as a data frame
  f <- read.table(i, sep = "\t", header = TRUE)

  # rename the columns
  colnames(f) <- c("gene_id", substr(i, 1, nchar(i)-7))

  #if the m is empty just copy the f to m
  if(length(m) == 0){
    m <- f
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    m <- merge(m, f, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f)
}

#grab the rows from the 1st column and use it as the row-names in the dataframe
rownames(m) <- m[,1]

# remove the column-1 (gene_ids) from the data frame using dplyr::select function
m <- m[,-1]

# the end result is a data frame with four columns containing the count data for each sample
# the row names are transcript IDs. 


#################################################
# Create metadata objects
#################################################

Sample = c("K32", "K23", "U13", "U32")
Condition = c("Killingworth", "Killingworth", "UConn", "UConn")
TimePoint = c("T2", "T3", "T3", "T2")
myfactors <- data.frame(Sample, Condition, TimePoint)
myfactors

# order and identity of samples in data frame m and myfactor must be the same

# ensure all the columns in dataframe m and myfactor is present
all(colnames(m) %in% myfactors$Sample)

# Then check the order is the same in both objects
all(colnames(m) == myfactors$Sample[1])

# They are not, so order them according to the sample names
m <- m[, myfactors$Sample]

# Now check the order is correct after sorting 
all(colnames(m) == myfactors$Sample)

# import transcript length information to a data frame
df_length <- read.table("length.txt", sep = "\t", header = TRUE, row.names = 1)

# create a vector to hold the length information
mylength <- df_length$length
names(mylength) <- row.names(df_length)
mylength <- mylength[rownames(m)]
head(mylength)


#################################################
# Preliminary data exploration and filtering
#################################################

# how many transcripts did we quantify expression for?
dim(m)
# we have 70k transcripts. 
# we tried to collapse them down to genes
# really 30k genes is an upper bound, so it's a bit unrealistic

# how many transcripts have zero expression for n samples?
zeroexpr <- (m==0) %>% rowSums()
table(zeroexpr)
# 24k transcripts expressed in only 1 sample!
# that's possibly concerning. 

# how much of the expression data maps to these transcripts?
colSums(m[zeroexpr==3,])/colSums(m)
# this is now quite problematic
# 30% of all data for k32 maps to transcripts unique to that sample compared to 1% or less for the other three. 
# this suggests a major issue with this sample. 
# we'd need to dig much further to figure out what that is. 


# what is the distribution of total expression across samples?
rowSums(m) %>% log(.,10) %>% hist(.,100)

# Before moving forward, we're going to exclude genes that are most likely to be useless in the analysis
# let's throw out all the transcripts with
  # a) less than 20 fragments mapping in total across all samples
  # b) that have non-zero expression for only a single sample

subg <- zeroexpr < 3 & rowSums(m) > 20
table(subg)

m <- m[subg,]
mylength <- mylength[subg]

#################################################
# Creating the NOISeq data object
#################################################

# Creating a NOISeq object
mydata <- readData(data = m, factors = myfactors, length = mylength)
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


# plot count distributions by sample:
countsbio <- dat(mydata, factor = NULL, type = "countsbio")
explo.plot(countsbio, toplot = 1, samples = NULL, plottype = "boxplot")

# do we need to apply among-sample normalization (the answer is almost always yes!)
mycd = dat(mydata, type = "cd", norm = FALSE, refColumn = 1)

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

mynoiseq.bio <- noiseqbio(mydata,
                          factor = "TimePoint",
                          norm = "tmm",
                          lc = 1)


res <- mynoiseq.bio@results[[1]]
head(res)

# MA plot

logmean <- log(rowSums(res[,1:2]),10)
logFC <- res$log2FC
sig <- res$prob >= 0.999
plot(logmean,logFC,pch=20,cex=.3,col=sig+1)

# let's make a heatmap of some DE genes:

# get log-transformed (+1) normalized counts
myTMM <- (tmm(assayData(mydata)$exprs, long = 1000, lc = 1)+1) %>% log(.,10)

# get 100 genes with prob(DE)=1
top100 <- which(sig)[1:100]
df <- myfactors[,-1]
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
colnames(mds)[4:5] <- c("mds1","mds2")

p <- ggplot(mds, aes(mds1, mds2, label = Sample)) +
  geom_point(aes(color = factor(TimePoint)))

p2 <- p + geom_text_repel(size=2) + labs(title = "mds plot of expression data")
p2

###############################################
# gene ontology enrichment analysis 
# using gfold top ranked genes
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



