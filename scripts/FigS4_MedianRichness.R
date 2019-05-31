# Teresita M. Porter, May 30, 2019

library(stringr)
library(scales)
library(reshape2)
library(vegan)
library(RColorBrewer)
library(ggplot2)
library("ggpubr")

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Filter table for Arthropoda only, work with ESVs only
B<-A[A$Phylum=="Arthropoda",]

#######################################################
# Create dataframe for rarefaction curves
## pool data across reps
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

# pivot to make esv matrix only
esv<-dcast(B3, kit+primer+site ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
esv$sample<-paste(esv$kit, esv$primer, esv$site, sep="_")

# delete duplicate columns
esv$kit<-NULL
esv$primer<-NULL
esv$site<-NULL

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

### Only do this for ESVs to work with taxa at the finest level of resolution ###
#Rarefy original ESV matrix down to 15th percentile library size to normalize read depth across samples

df<-rrarefy(esv_notnull2, sample=esv_15percentile)

###################################################################
##### Calculate richness
###################################################################

# Get total number of species per sample
# ## Separate out soil and tissue DNA extraction kit results

## For each kit, calculate specnumber (number of ESVs) for each primer
df_soil_A<-specnumber(data.frame(df[grep("^SoilKit_A_",rownames(df)),]))
df_soil_B<-specnumber(data.frame(df[grep("^SoilKit_B_",rownames(df)),]))
df_soil_C<-specnumber(data.frame(df[grep("^SoilKit_C_",rownames(df)),]))
df_soil_D<-specnumber(data.frame(df[grep("^SoilKit_D_",rownames(df)),]))
df_soil_E<-specnumber(data.frame(df[grep("^SoilKit_E_",rownames(df)),]))
df_soil_F<-specnumber(data.frame(df[grep("^SoilKit_F_",rownames(df)),]))
#df_tissue_A<-specnumber(data.frame(df[grep("^TissueKit_A_",rownames(df)),]))
#df_tissue_B<-specnumber(data.frame(df[grep("^TissueKit_B_",rownames(df)),]))

# Create matrix for richness for primer comparison, do this for each kit
richness<-data.frame(c(df_soil_A, df_soil_B, df_soil_C, df_soil_D, df_soil_E, df_soil_F))
names(richness)[1]<-"richness"

# Create df for ggplot using rownames for soil primer expt
## Grab row names and make new column
richness$samples<-row.names(richness)

## Split samples column into their own fields
richness[,3:5]<-do.call('rbind', strsplit(as.character(richness$samples),'_',fixed=TRUE))
## Remove samples column
richness<-richness[,-2]
## Rename columns
names(richness)[2:4]<-c("kit","primer","site")

#create factors and levels for soil and tissue richness df's
richness$primer<-factor(richness$primer, 
                        levels=c("A","B","C","D","E","F"),
                        labels=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"))
richness$kit<-factor(richness$kit,
                     levels=c("SoilKit"),
                     labels=c("Soil Kit"))

# Create custom color palette, remove yellow because hard to see
custom<-brewer.pal(n = 7, name = 'Set1')
custom2<-custom[-6]

#Visualize with ggplot
p<-ggplot(richness, aes(x=primer,y=richness)) +
  geom_boxplot(alpha=0.8) +
  geom_point(aes(color=primer), size=2, shape=20, position=position_jitterdodge()) +
  labs(x="Primer", y="ESV Richness") +
  theme(text = element_text(size=16),
        axis.text.x = element_text(angle=90, hjust=1),
        plot.title=element_text(hjust=0),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key=element_blank(),
        legend.title=element_blank()) +
  guides(color=FALSE) +
  scale_color_manual(values=custom2)
#  facet_wrap(~ kit, scales="free_x")

ggsave("FigS4_MedianRichness.pdf",p)

###################################
# Check for normality
###################################

# Visual tests for normality
ggdensity(richness$richness, 
          main = "density richness dist",
          xlab = "richness")

ggqqplot(richness$richness)

# Shapiro-Wilk test of normality (often positive for small sample sizes)
shapiro.test(richness$richness)
# W = 0.97388, p-value = 0.3558 NOT SIGNIFICANTLY different from normal distributions

# pairwise t-tests with Holm adjustment for multiple comparisons for SOIL kit
pairwise.t.test(richness$richness[richness$kit=="Soil Kit"], richness$primer[richness$kit=="Soil Kit"], p.adj="holm")
# BR5   F230R ml-jg BF1R2 BF2R2
# F230R 1.000 -     -     -     -    
#   ml-jg 1.000 1.000 -     -     -    
#   BF1R2 1.000 1.000 1.000 -     -    
#   BF2R2 1.000 1.000 0.461 1.000 -    
#   fwh1  0.299 0.209 0.059 0.299 1.000
