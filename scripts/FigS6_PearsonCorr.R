# Teresita M. Porter, May 31, 2019

library(stringr)
library(scales)
library(reshape2)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library("ggpubr")
library(corrplot)
library(data.table)
library(gridExtra)
library(grid)
library("psych") #corr.test # remove

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Filter table for Arthropoda only, work with ESVs only
B<-A[A$Phylum=="Arthropoda",]

#######################################################
# Create dataframe for rarefaction curves
## DO NOT pool data across reps
## separate curves for primers and kits and sites
######################################################

# Split up SampleName for Arthropoda matrix (for esv, order, class ranks)
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
# Add new column names
names(B2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
# create df for each primer and add primer name to a new primer column
B2_A<-B2[grep("^A",B2$A_GlobalESV),]
B2_A$primer<-"A"
B2_B<-B2[grep("^B",B2$A_GlobalESV),]
B2_B$primer<-"B"
B2_C<-B2[grep("^C",B2$A_GlobalESV),]
B2_C$primer<-"C"
B2_D<-B2[grep("^D",B2$A_GlobalESV),]
B2_D$primer<-"D"
B2_E<-B2[grep("^E",B2$A_GlobalESV),]
B2_E$primer<-"E"
B2_F<-B2[grep("^F",B2$A_GlobalESV),]
B2_F$primer<-"F"
#combine them all into a single dataframe
B3<-rbind(B2_A, B2_B, B2_C, B2_D, B2_E, B2_F)

# pivot to make esv matrix only (keep reps separate here)
esv<-dcast(B3, kit+primer+site+rep ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
esv$sample<-paste(esv$kit, esv$primer, esv$site, esv$rep, sep="_")

# delete duplicate columns
esv$kit<-NULL
esv$primer<-NULL
esv$site<-NULL
esv$rep<-NULL

#move SampleName to rownames then delete
rownames(esv)<-esv$sample
esv$sample<-NULL

#remove columns with only zeros
esv_notnull<-esv[,colSums(esv) !=0]

#remove rows with only zeros & edit rownames
esv_notnull2<-esv_notnull[rowSums(esv_notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv_15percentile<-quantile(rowSums(esv_notnull2), prob=0.15)

###################################################################
##### Rarefy dataset
###################################################################

set.seed(1234)

df<-rrarefy(esv_notnull2, sample=esv_15percentile)

##############################################
#### Convert to presence-absence matrix
##############################################

df[df>0]<-1

# rep test 
c<-df[(grepl("^Soil", rownames(df))) &
      (grepl("_1$", rownames(df))),]
d<-df[(grepl("^Soil", rownames(df))) &
      ( grepl("_2$", rownames(df))),]

# transpose to get interesting quantities (samples in columns)
c_t<-t(c)
d_t<-t(d)

##############################################
#### Create correlation plot
##############################################

pdf("FigS6_corplot_rep.pdf")

# calc correlations bewtween PCR replicates
rep_pearson<-corr.test(c_t,d_t, method="pearson", adjust="holm", use="pairwise")

corrplot(rep_pearson$r, p.mat = rep_pearson$p, sig.level=0.05, insig="blank",
        method="circle", order="hclust", hclust.method="average", type="upper", 
        col = brewer.pal(n = 8, name = "RdBu"), tl.col="black", cl.ratio = 0.5, cl.align = "r", 
        cl.cex=1, tl.cex=0.7)

dev.off()
