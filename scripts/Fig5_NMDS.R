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
library(plyr)
library(pairwiseAdonis)

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Filter table for Arthropoda only, work with ESVs only
B<-A[A$Phylum=="Arthropoda",]

# NMDS at ESV rank just clusters by primer, try species rank sBP>=0.70 (95% correct)
C<-B[B$sBP>=0.70,]

# reassign to B to make code changes easier when trying different ranks
B<-C

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

# Only work with SoilKit (ignore TissueKit)
B4<-B3[B3$kit=="SoilKit",]

# pivot to make esv matrix only (keep reps separate here)
esv<-dcast(B4, primer+site+rep ~ Species, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
esv$sample<-paste(esv$primer, esv$site, esv$rep, sep="_")

# delete duplicate columns
#esv$kit<-NULL
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

#Rarefy original ESV matrix down to 15th percentile library size to normalize read depth across samples

df<-rrarefy(esv_notnull2, sample=esv_15percentile)

##############################################
#### Convert to presence-absence matrix
##############################################

df[df>0]<-1

###################################################################
##### Create distance matrix for beta diversity analysis at ESV rank
###################################################################

#Do 2 dimensional NMDS, no environment file for now 
nmds2<-metaMDS(df, k=2,trymax=100)
# k=2, stress=0.154

# Create grouping matrix for samples
## Grab row names from matrix above
sample_df<-data.frame(row.names(df))
## Rename the column
names(sample_df)<-"rn"
## Copy column to row names
row.names(sample_df)<-sample_df$rn
## Split first column into their own fields
sample_df[,2:4]<-do.call('rbind', strsplit(as.character(sample_df$rn),'_',fixed=TRUE))
## Remove first column
sample_df<-sample_df[,-1]
## Rename columns
names(sample_df)<-c("primer","site","rep")

# Grab sites & scores
#A_spp.sc <- scores(ILC_nmds2, display = "species")
site.sc <- data.frame(scores(nmds2, display = "sites"))

# Put it all in one df
merged <- merge(site.sc,sample_df,by="row.names")
colnames(merged)[colnames(merged)=="Row.names"] <- "rn"

#create factors and levels for chulls
merged$site<-factor(merged$site, levels=c("1","2","3","4","5","6"))
merged$primer<-factor(merged$primer, levels=c("A","B","C","D","E","F"),
                      labels=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"))

# Create convex hulls
chulls <- ddply(merged, .(site), function(merged) merged[chull(merged$NMDS1, merged$NMDS2), ])

# list of markers
markers=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1")

# Create scatter plot
p1<-ggplot(data=merged, aes(x=NMDS1, y=NMDS2)) + 
  geom_polygon(data=chulls, aes(x=NMDS1, y=NMDS2, fill=site), alpha=1) +
  geom_text(data=merged, aes(x=NMDS1, y=NMDS2, label=primer), size=3) +
  scale_fill_manual(name="Site",values=c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")) +
  theme_bw() +
  theme(text = element_text(size=16),
        axis.text.x = element_text(hjust=1),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_blank()) +
  guides(fill = guide_legend(order = 1), 
         colour = FALSE, 
         shape = FALSE)

ggsave("Fig5_NMDS.pdf", p1)

##############################
#Assess dispersion (variance) using ANOVA
##############################

# Create distance matrix based on P-A data using Bray/Sorensen dissimilarity
sor<-vegdist(df, "bray", binary=TRUE)

# Calculate beta dispersion (homogeneity needed for adonis)
#primer
bd_primer<-betadisper(sor, factor(merged[,4]))
#site
bd_site<-betadisper(sor, factor(merged[,5]))
#rep
bd_rep<-betadisper(sor, factor(merged[,6]))

# check for homogeneity of beta dispersions within groups
anova(bd_primer) # n/s
anova(bd_site) # 0.03776 * heterogenous beta disp, but we have a balanced design
anova(bd_rep) # n/s

pdf("BetaDispersion_boxplot.pdf")
par(mfrow=c(2,2))
boxplot(bd_primer, las=2)
mtext("Primers", side=1, line=4 )
mtext("A)", side=3, adj=0, line=1.2)
boxplot(bd_site)
mtext("Sites", side=1, line=4)
mtext("B)", side=3, adj=0, line=1.2)
boxplot(bd_rep)
mtext("Replicates", side=1, line=4)
mtext("C)", side=3, adj=0, line=1.2)
dev.off()

###################################################################
##### Shephards curve and goodness of fit calcs
###################################################################

# Stressplot Shephards curve to assess goodness of fit between observed and ordination distances
pdf("stressplot.pdf") # Linear fit, R2 = 0.912
stressplot(nmds2)
gof <-goodness(nmds2)
gof
plot(nmds2, display = "sites", type="n")
points(nmds2, display="sites",cex=2*gof/mean(gof))
dev.off()

############
# Use ADONIS to test for significant interactions between groups
############

adonis(df~site*primer*rep, data=merged, permutations=999)
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# site             5   14.5141       3       0 0.75991      1
# primer           5    1.3308       0       0 0.06967      1
# rep              1    0.0146       0       0 0.00076      1
# site:primer     25    2.1685       0       0 0.11354      1
# site:rep         5    0.0718       0       0 0.00376      1
# primer:rep       5    0.1018       0       0 0.00533      1
# site:primer:rep 25    0.8983       0       0 0.04703      1
# Residuals        0    0.0000     Inf         0.00000       
# Total           71   19.0999                 1.00000       

# no significant interactinos found
adonis(df~site, data=merged, permutations=999)
# 0.001 *** sites are significantly diff; R2=0.75959 where 76% variation explained by site

adonis(df~primer, data=merged, permutations=999)
# n/s

adonis(df~rep, data=merged, permutations=999)
# n/s
