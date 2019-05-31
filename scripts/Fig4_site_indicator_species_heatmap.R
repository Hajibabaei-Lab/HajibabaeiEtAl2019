# Teresita M. Porter, May 30, 2019

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
library(indicspecies)
library(janitor)

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Filter table for Arthropoda only, use for ESVs, order, class
B<-A[A$Phylum=="Arthropoda",]

# For species use sBP>=0.70 (95% correct)
C<-B[B$sBP>=0.70,]

#######################################################
# Create dataframe for SPECIES rarefaction
## DO NOT pool data across reps
######################################################

# Split up SampleName for Arthropoda matrix, filtered by sBP>=0.70
C2<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")))

# Add new column names
names(C2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
C2_A<-C2[grep("^A",C2$A_GlobalESV),]
C2_A$primer<-"A"
C2_B<-C2[grep("^B",C2$A_GlobalESV),]
C2_B$primer<-"B"
C2_C<-C2[grep("^C",C2$A_GlobalESV),]
C2_C$primer<-"C"
C2_D<-C2[grep("^D",C2$A_GlobalESV),]
C2_D$primer<-"D"
C2_E<-C2[grep("^E",C2$A_GlobalESV),]
C2_E$primer<-"E"
C2_F<-C2[grep("^F",C2$A_GlobalESV),]
C2_F$primer<-"F"

#combine them all into a single dataframe
C3<-rbind(C2_A, C2_B, C2_C, C2_D, C2_E, C2_F)

# pivot to make SPECIES matrix only (keep reps separate here)
## proceed processing with species (because multipatt doesn't recover as many indicators when the extra long names are used, truncated during processing?)
species<-dcast(C3, kit+primer+site+rep ~ Species, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
species$sample<-paste(species$kit, species$primer, species$site, species$rep, sep="_")

# delete duplicate columns
species$kit<-NULL
species$primer<-NULL
species$site<-NULL
species$rep<-NULL

#move SampleName to rownames then delete
rownames(species)<-species$sample
species$sample<-NULL

#remove columns with only zeros
species_notnull<-species[,colSums(species) !=0]

#remove rows with only zeros & edit rownames
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]

#calculate 15th percentile for rrarefy function
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)

# Rarefy dataset
set.seed(1234)
df<-rrarefy(species_notnull2, sample=species_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Do indicator SPECIES analysis
# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset  the markers, work with soil only, pool both replicates
Aprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_A_", rownames(df.2)),]
Bprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_B_", rownames(df.2)),]
Cprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_C_", rownames(df.2)),]
Dprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_D_", rownames(df.2)),]
Eprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_E_", rownames(df.2)),]
Fprimer<-df.2[grepl("^Soil", rownames(df.2)) &
          grepl("_F_", rownames(df.2)),]

# Remove any zero-sum columns
A_notnull<-Aprimer[,colSums(Aprimer) !=0]
B_notnull<-Bprimer[,colSums(Bprimer) !=0]
C_notnull<-Cprimer[,colSums(Cprimer) !=0]
D_notnull<-Dprimer[,colSums(Dprimer) !=0]
E_notnull<-Eprimer[,colSums(Eprimer) !=0]
F_notnull<-Fprimer[,colSums(Fprimer) !=0]

# Sort out groups (sites) for each marker 
## do rownames(A_notnull), groups same for each marker (balanced design)
groups<-c(1,1,2,2,3,3,4,4,5,5,6,6)

# Do indicator species analysis
A_indval=multipatt(A_notnull, groups, control=how(nperm=999))
B_indval=multipatt(B_notnull, groups, control=how(nperm=999))
C_indval=multipatt(C_notnull, groups, control=how(nperm=999))
D_indval=multipatt(D_notnull, groups, control=how(nperm=999))
E_indval=multipatt(E_notnull, groups, control=how(nperm=999))
F_indval=multipatt(F_notnull, groups, control=how(nperm=999))

# Extract species list with pvalues from each analysis
A_indval_df<-data.frame(A_indval$sign)
B_indval_df<-data.frame(B_indval$sign)
C_indval_df<-data.frame(C_indval$sign)
D_indval_df<-data.frame(D_indval$sign)
E_indval_df<-data.frame(E_indval$sign)
F_indval_df<-data.frame(F_indval$sign)

# Grab indicator species with pvalues <=0.05
A_indic<-A_indval_df[which(A_indval_df$p.value<=0.05),]
B_indic<-B_indval_df[which(B_indval_df$p.value<=0.05),]
C_indic<-C_indval_df[which(C_indval_df$p.value<=0.05),]
D_indic<-D_indval_df[which(D_indval_df$p.value<=0.05),]
E_indic<-E_indval_df[which(E_indval_df$p.value<=0.05),]
F_indic<-F_indval_df[which(F_indval_df$p.value<=0.05),]

# Move row names to first column 'rn' using data.table
setDT(A_indic, keep.rownames = TRUE)[]
setDT(B_indic, keep.rownames = TRUE)[]
setDT(C_indic, keep.rownames = TRUE)[]
setDT(D_indic, keep.rownames = TRUE)[]
setDT(E_indic, keep.rownames = TRUE)[]
setDT(F_indic, keep.rownames = TRUE)[]

# Pool data across sites
A_indic$pooled<-A_indic$s.1 + A_indic$s.2 + A_indic$s.3 + A_indic$s.4 + A_indic$s.5 + A_indic$s.6
B_indic$pooled<-B_indic$s.1 + B_indic$s.2 + B_indic$s.3 + B_indic$s.4 + B_indic$s.5 + B_indic$s.6
C_indic$pooled<-C_indic$s.1 + C_indic$s.2 + C_indic$s.3 + C_indic$s.4 + C_indic$s.5 + C_indic$s.6
D_indic$pooled<-D_indic$s.1 + D_indic$s.2 + D_indic$s.3 + D_indic$s.4 + D_indic$s.5 + D_indic$s.6
E_indic$pooled<-E_indic$s.1 + E_indic$s.2 + E_indic$s.3 + E_indic$s.4 + E_indic$s.5 + E_indic$s.6
F_indic$pooled<-F_indic$s.1 + F_indic$s.2 + F_indic$s.3 + F_indic$s.4 + F_indic$s.5 + F_indic$s.6

# Convert all rn to presence absence (except for Total)
A_indic$pooled[(A_indic$pooled>0)]<-1
B_indic$pooled[(B_indic$pooled>0)]<-1
C_indic$pooled[(C_indic$pooled>0)]<-1
D_indic$pooled[(D_indic$pooled>0)]<-1
E_indic$pooled[(E_indic$pooled>0)]<-1
F_indic$pooled[(F_indic$pooled>0)]<-1

# Add a row of totals
A_indic<-adorn_totals(A_indic, where="row")
B_indic<-adorn_totals(B_indic, where="row")
C_indic<-adorn_totals(C_indic, where="row")
D_indic<-adorn_totals(D_indic, where="row")
E_indic<-adorn_totals(E_indic, where="row")
F_indic<-adorn_totals(F_indic, where="row")

# Add column to indicate marker
A_indic$primer<-"A"
B_indic$primer<-"B"
C_indic$primer<-"C"
D_indic$primer<-"D"
E_indic$primer<-"E"
F_indic$primer<-"F"

# Consolidate these results into single df for plotting
merged<- rbind(A_indic, B_indic, C_indic, D_indic, E_indic, F_indic)

# Change NA to 0
merged[is.na(merged)] <- 0

# Change 'Total' in merged$rn to 'ZTotal'
merged$rn[merged$rn=="Total"]<-"ZTotal"

# Create factors
merged$primer<-factor(merged$primer, levels=c("A","B","C","D","E","F"),
                      labels=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"))

# copy df then edit names
merged2<-merged

# Create mapping file
# C3 is the taxonomy file filtered for sBP>=0.70
# Just keep the Order, Family, and Species fields, uniques only
map<-unique(C3[,c("Order", "Family", "Species")])
# combine Order, Family, and Species into a single lineage field
map$lineage<-paste(map$Order, map$Family, map$Species)
# remove unneeded Order and Family fields
# Keep Species and lineage only as the mapping file
map<-map[,-c(1:2)]
# Account for Ztotal
new<-data.frame("Species"="ZTotal","lineage"="ZTotal")
map<-rbind(map, new)

# create function to replace Species with lineage
merged2$rn <- map$lineage[match(unlist(merged$rn), map$Species)]

# Change underscores to space
merged2$rn<-gsub("_", " ", merged2$rn)

# Create heat map
labs <- sapply(strsplit(as.character(merged2$rn), " "), 
               FUN = function(x) {
                 x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4];
                 parse(text = paste("plain('", x1, "') ~ ", "plain('", x2, "') ~ ", "italic('", x3, "') ~ ", "italic('", x4, "')", sep = ""))
               })

p2<-ggplot(data=merged2, aes(x=primer,y=rn)) +
  geom_tile(data=merged2[merged2$pooled>=1]) +
  coord_equal() +
  scale_y_discrete(limits = rev(sort(unique(merged2$rn))),
                   labels=labs,
                   breaks=merged2$rn) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_text(size=7),
        axis.text.x=element_text(size=7, angle=90),
        legend.text=element_text(size=7),
        legend.title=element_text(size=7),
        legend.key.size = unit(.5,"line")) +
  annotate("text", x=1, y=1, label="14", size=2.5) +
  annotate("text", x=2, y=1, label="11", size=2.5) +
  annotate("text", x=3, y=1, label="10", size=2.5) +
  annotate("text", x=4, y=1, label="11", size=2.5) +
  annotate("text", x=5, y=1, label="13", size=2.5) +
  annotate("text", x=6, y=1, label="5", size=2.5)

ggsave("F4_site_indicator_species_heatmap.pdf")
