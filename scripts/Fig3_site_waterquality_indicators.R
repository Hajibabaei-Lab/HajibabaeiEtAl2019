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

# Filter table for Arthropoda only, use for arthropoda site indicators at ESVs
B<-A[A$Phylum=="Arthropoda",]

# For species use sBP>=0.70 (95% correct), use for arthropoda site indicators at species rank
C<-B[B$sBP>=0.70,]

# For genera use gBP>=0.30 (99% correct), use for arthropoda site indicators at genus rank
D<-B[B$gBP>=0.30,]

# For family use fBP>=0.20 (99% correct), use for arthropoda site indicators and EPTC indic at family rank
E<-B[B$fBP>=0.20,]

# Use with EPTC at ESV rank
F<-data.frame(matrix(ncol=ncol(B), nrow=0))
names(F) <- colnames(B)

# exclude Chironomidae with low support, keep the rest
for (i in 1:nrow(B)) {
  if (B[i,23]=="Chironomidae") { #Family
    if (B[i,25]>=0.20) { #fBP
      F<-rbind(F,B[i,])
    } 
  } 
  else {
    F <- rbind(F,B[i,])
  }
}

# Use with EPTC at species rank
G<-data.frame(matrix(ncol=ncol(B), nrow=0))
names(G) <- colnames(B)

# exclude Chironomidae with low support, exclude species with low support
for (i in 1:nrow(B)) {
  if (B[i,23]=="Chironomidae") { #Family
    if (B[i,25]>=0.20) { #fBP
      if (B[i,31]>=0.70) { #sBP
        G<-rbind(G,B[i,])
      }
    } 
  } 
  else {
    if (B[i,31]>=0.70) { #sBP
      G <- rbind(G,B[i,])      
    }
  }
}

# Use with EPTC at genus rank
H<-data.frame(matrix(ncol=ncol(B), nrow=0))
names(H) <- colnames(B)

# exclude Chironomidae with low support, exclude species with low support
for (i in 1:nrow(B)) {
  if (B[i,23]=="Chironomidae") { #Family
    if (B[i,25]>=0.20) { #fBP
      if (B[i,31]>=0.30) { #gBP
        H<-rbind(H,B[i,])
      }
    } 
  } 
  else {
    if (B[i,31]>=0.30) { #gBP
      H <- rbind(H,B[i,])      
    }
  }
}

###############################################################################
# Arthropoda site indicator analysis
###############################################################################

#######################################################
# Create dataframe for ESV rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
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

# pivot to make ESV matrix only
esv<-dcast(B3, kit+primer+site+rep ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+primer into single column
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

# Rarefy dataset
set.seed(1234)
df<-rrarefy(esv_notnull2, sample=esv_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset the markers, work with soil only, pool both replicates
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
groups<-c(1,1,2,2,3,3,4,4,5,5,6,6)

# Do Arhtropoda site indicator ESV analysis
A_indval=multipatt(A_notnull, groups, control=how(nperm=999))
B_indval=multipatt(B_notnull, groups, control=how(nperm=999))
C_indval=multipatt(C_notnull, groups, control=how(nperm=999))
D_indval=multipatt(D_notnull, groups, control=how(nperm=999))
E_indval=multipatt(E_notnull, groups, control=how(nperm=999))
F_indval=multipatt(F_notnull, groups, control=how(nperm=999))

# Extract ESVs list with pvalues from each analysis
A_indval_df<-data.frame(A_indval$sign)
B_indval_df<-data.frame(B_indval$sign)
C_indval_df<-data.frame(C_indval$sign)
D_indval_df<-data.frame(D_indval$sign)
E_indval_df<-data.frame(E_indval$sign)
F_indval_df<-data.frame(F_indval$sign)

# Grab indicator ESVs with pvalues <=0.05
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

# Convert all rn to presence absence
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

# Just create a plot that summarizes the results

esv_summary<-merged[merged$rn=="Total"]
esv_summary$rank<-"esv"

#######################################################
# Create dataframe for SPECIES rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
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

# pivot to make SPECIES matrix only 
species<-dcast(C3, kit+primer+site+rep ~ Species, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+primer+site+rep into single column
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

#remove rows with only zeros
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]

#calculate 15th percentile for rrarefy function
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)

# Rarefy dataset
set.seed(1234)
df<-rrarefy(species_notnull2, sample=species_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

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

# Do Arthropoda site indicator species analysis
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

# Convert all rn to presence absence
A_indic$pooled[(A_indic$pooled>0)]<-1
B_indic$pooled[(B_indic$pooled>0)]<-1
C_indic$pooled[(C_indic$pooled>0)]<-1
D_indic$pooled[(D_indic$pooled>0)]<-1
E_indic$pooled[(E_indic$pooled>0)]<-1
F_indic$pooled[(F_indic$pooled>0)]<-1

# Add column totals
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

# Create a plot that summarizes the results
species_summary<-merged[merged$rn=="Total"]
species_summary$rank<-"species"

#######################################################
# Create dataframe for Genus rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
D2<-data.frame(D, do.call(rbind, str_split(D$SampleName,"_")))

# Add new column names
names(D2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
D2_A<-D2[grep("^A",D2$A_GlobalESV),]
D2_A$primer<-"A"
D2_B<-D2[grep("^B",D2$A_GlobalESV),]
D2_B$primer<-"B"
D2_C<-D2[grep("^C",D2$A_GlobalESV),]
D2_C$primer<-"C"
D2_D<-D2[grep("^D",D2$A_GlobalESV),]
D2_D$primer<-"D"
D2_E<-D2[grep("^E",D2$A_GlobalESV),]
D2_E$primer<-"E"
D2_F<-D2[grep("^F",D2$A_GlobalESV),]
D2_F$primer<-"F"

#combine them all into a single dataframe
D3<-rbind(D2_A, D2_B, D2_C, D2_D, D2_E, D2_F)

# pivot to make genus matrix only 
genus<-dcast(D3, kit+primer+site+rep ~ Genus, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
genus$sample<-paste(genus$kit, genus$primer, genus$site, genus$rep, sep="_")

# delete duplicate columns
genus$kit<-NULL
genus$primer<-NULL
genus$site<-NULL
genus$rep<-NULL

#move SampleName to rownames then delete
rownames(genus)<-genus$sample
genus$sample<-NULL

#remove columns with only zeros
genus_notnull<-genus[,colSums(genus) !=0]

#remove rows with only zeros & edit rownames
genus_notnull2<-genus_notnull[rowSums(genus_notnull) !=0,]

#calculate 15th percentile for rrarefy function
genus_15percentile<-quantile(rowSums(genus_notnull2), prob=0.15)

# Rarefy dataset
set.seed(1234)
df<-rrarefy(genus_notnull2, sample=genus_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

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

# Do Arthropoda site indicator genus analysis
A_indval=multipatt(A_notnull, groups, control=how(nperm=999))
B_indval=multipatt(B_notnull, groups, control=how(nperm=999))
C_indval=multipatt(C_notnull, groups, control=how(nperm=999))
D_indval=multipatt(D_notnull, groups, control=how(nperm=999))
E_indval=multipatt(E_notnull, groups, control=how(nperm=999))
F_indval=multipatt(F_notnull, groups, control=how(nperm=999))

# Extract genuss list with pvalues from each analysis
A_indval_df<-data.frame(A_indval$sign)
B_indval_df<-data.frame(B_indval$sign)
C_indval_df<-data.frame(C_indval$sign)
D_indval_df<-data.frame(D_indval$sign)
E_indval_df<-data.frame(E_indval$sign)
F_indval_df<-data.frame(F_indval$sign)

# Grab indicator genuss with pvalues <=0.05
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

# Convert all rn to presence absence
A_indic$pooled[(A_indic$pooled>0)]<-1
B_indic$pooled[(B_indic$pooled>0)]<-1
C_indic$pooled[(C_indic$pooled>0)]<-1
D_indic$pooled[(D_indic$pooled>0)]<-1
E_indic$pooled[(E_indic$pooled>0)]<-1
F_indic$pooled[(F_indic$pooled>0)]<-1

# Add row of totals
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

# Just create a plot that summarizes the results
genus_summary<-merged[merged$rn=="Total"]
genus_summary$rank<-"genus"

#######################################################
# Create dataframe for Family rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
E2<-data.frame(E, do.call(rbind, str_split(E$SampleName,"_")))

# Add new column names
names(E2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
E2_A<-E2[grep("^A",E2$A_GlobalESV),]
E2_A$primer<-"A"
E2_B<-E2[grep("^B",E2$A_GlobalESV),]
E2_B$primer<-"B"
E2_C<-E2[grep("^C",E2$A_GlobalESV),]
E2_C$primer<-"C"
E2_D<-E2[grep("^D",E2$A_GlobalESV),]
E2_D$primer<-"D"
E2_E<-E2[grep("^E",E2$A_GlobalESV),]
E2_E$primer<-"E"
E2_F<-E2[grep("^F",E2$A_GlobalESV),]
E2_F$primer<-"F"

#combine them all into a single dataframe
E3<-rbind(E2_A, E2_B, E2_C, E2_D, E2_E, E2_F)

# pivot to make family matrix only (keep reps separate here)
family<-dcast(E3, kit+primer+site+rep ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
family$sample<-paste(family$kit, family$primer, family$site, family$rep, sep="_")

# delete duplicate columns
family$kit<-NULL
family$primer<-NULL
family$site<-NULL
family$rep<-NULL

#move SampleName to rownames then delete
rownames(family)<-family$sample
family$sample<-NULL

#remove columns with only zeros
family_notnull<-family[,colSums(family) !=0]

#remove rows with only zeros & edit rownames
family_notnull2<-family_notnull[rowSums(family_notnull) !=0,]

#calculate 15th percentile for rrarefy function
family_15percentile<-quantile(rowSums(family_notnull2), prob=0.15)

# Rarefy dataset
set.seed(1234)
df<-rrarefy(family_notnull2, sample=family_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

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

# Do Arthropoda site indicator family analysis
A_indval=multipatt(A_notnull, groups, control=how(nperm=999))
B_indval=multipatt(B_notnull, groups, control=how(nperm=999))
C_indval=multipatt(C_notnull, groups, control=how(nperm=999))
D_indval=multipatt(D_notnull, groups, control=how(nperm=999))
E_indval=multipatt(E_notnull, groups, control=how(nperm=999))
F_indval=multipatt(F_notnull, groups, control=how(nperm=999))

# Extract familys list with pvalues from each analysis
A_indval_df<-data.frame(A_indval$sign)
B_indval_df<-data.frame(B_indval$sign)
C_indval_df<-data.frame(C_indval$sign)
D_indval_df<-data.frame(D_indval$sign)
E_indval_df<-data.frame(E_indval$sign)
F_indval_df<-data.frame(F_indval$sign)

# Grab indicator families with pvalues <=0.05
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

# Convert all rn to presence absence
A_indic$pooled[(A_indic$pooled>0)]<-1
B_indic$pooled[(B_indic$pooled>0)]<-1
C_indic$pooled[(C_indic$pooled>0)]<-1
D_indic$pooled[(D_indic$pooled>0)]<-1
E_indic$pooled[(E_indic$pooled>0)]<-1
F_indic$pooled[(F_indic$pooled>0)]<-1

# Add column totals
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

# Just create a plot that summarizes the results
family_summary<-merged[merged$rn=="Total"]
family_summary$rank<-"family"

# Put all the summry results in ONE table
summary<-rbind(esv_summary, species_summary, genus_summary, family_summary)

# Add column to indicate indicator type
summary$ind<-"broad"

# remove unneeded columns
summary<-summary[,-c(1,2,3,4,5,6,7,8,9,10)]

###############################################################################
# Freshwater indicator analysis
###############################################################################

#######################################################
# Create dataframe for ESV rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
F2<-data.frame(F, do.call(rbind, str_split(F$SampleName,"_")))

# Add new column names
names(F2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
F2_A<-F2[grep("^A",F2$A_GlobalESV),]
F2_A$primer<-"A"
F2_B<-F2[grep("^B",F2$A_GlobalESV),]
F2_B$primer<-"B"
F2_C<-F2[grep("^C",F2$A_GlobalESV),]
F2_C$primer<-"C"
F2_D<-F2[grep("^D",F2$A_GlobalESV),]
F2_D$primer<-"D"
F2_E<-F2[grep("^E",F2$A_GlobalESV),]
F2_E$primer<-"E"
F2_F<-F2[grep("^F",F2$A_GlobalESV),]
F2_F$primer<-"F"

#combine them all into a single dataframe
F3<-rbind(F2_A, F2_B, F2_C, F2_D, F2_E, F2_F)

# pivot to make ESV matrix only (pool data across sites + reps)
esv<-dcast(F3, kit+primer ~ A_GlobalESV+Order+Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
esv$sample<-paste(esv$kit, esv$primer, sep="_")

# delete duplicate columns
esv$kit<-NULL
esv$primer<-NULL

#move SampleName to rownames then delete
rownames(esv)<-esv$sample
esv$sample<-NULL

#remove columns with only zeros
esv_notnull<-esv[,colSums(esv) !=0]

#remove rows with only zeros & edit rownames
esv_notnull2<-esv_notnull[rowSums(esv_notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv_15percentile<-quantile(rowSums(esv_notnull2), prob=0.15)

# Normalize dataset down to 15th percentile
set.seed(1234)
df<-rrarefy(esv_notnull2, sample=esv_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset the markers, work with soil only, pool both replicates
df.2<-df.2[grepl("^Soil", rownames(df.2)),]

# Remove zero sum columns
df.2_notnull<-df.2[,colSums(df.2) !=0]

# Retain only the freshwater groups
df.2_notnull<-df.2_notnull[,grepl("_Ephemeroptera_|_Plecoptera_Insecta_|_Trichoptera_|_Chironomidae", names(df.2_notnull))]

# Trim colnames to just the ESV // keep first two fields, not remove last 2
names <- gsub("^([A-F]+_Otu\\d+)_\\w+_\\w+_*\\w*", "\\1", colnames(df.2_notnull))
names(df.2_notnull)<-names

# move rownames to first col
setDT(df.2_notnull, keep.rownames = TRUE)[]

# put primer into own column
df.2_notnull$primer <- str_sub(df.2_notnull$rn, -1,-1)

# add a column of totals
esv.df<-adorn_totals(df.2_notnull, where="col")

# put into summary table for ggplot
esv.df2<-data.frame(esv.df$Total, esv.df$primer)

#rename cols
names(esv.df2) <- c("pooled","primer")

# add column for rank
esv.df2$rank<-"esv"

# add column for indicator
esv.df2$ind<-"fresh"

############################################
# Create dataframe for Species rarefaction
############################################

# Split up SampleName for Arthropoda matrix
G2<-data.frame(G, do.call(rbind, str_split(G$SampleName,"_")))

# Add new column names
names(G2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
G2_A<-G2[grep("^A",G2$A_GlobalESV),]
G2_A$primer<-"A"
G2_B<-G2[grep("^B",G2$A_GlobalESV),]
G2_B$primer<-"B"
G2_C<-G2[grep("^C",G2$A_GlobalESV),]
G2_C$primer<-"C"
G2_D<-G2[grep("^D",G2$A_GlobalESV),]
G2_D$primer<-"D"
G2_E<-G2[grep("^E",G2$A_GlobalESV),]
G2_E$primer<-"E"
G2_F<-G2[grep("^F",G2$A_GlobalESV),]
G2_F$primer<-"F"

#combine them all into a single dataframe
G3<-rbind(G2_A, G2_B, G2_C, G2_D, G2_E, G2_F)

# pivot to make species matrix only (pool data across sites and reps)
species<-dcast(G3, kit+primer ~ Species+Order+Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
species$sample<-paste(species$kit, species$primer, sep="_")

# delete duplicate columns
species$kit<-NULL
species$primer<-NULL

#move SampleName to rownames then delete
rownames(species)<-species$sample
species$sample<-NULL

#remove columns with only zeros
species_notnull<-species[,colSums(species) !=0]

#remove rows with only zeros & edit rownames
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]

#calculate 15th percentile for rrarefy function
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)

# Normalize dataset down to 15th percentile
set.seed(1234)
df<-rrarefy(species_notnull2, sample=species_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset the markers, work with soil only
df.2<-df.2[grepl("^Soil", rownames(df.2)),]

# Remove zero sum columns
df.2_notnull<-df.2[,colSums(df.2) !=0]

# Retain only the freshwater groups
df.2_notnull<-df.2_notnull[,grepl("_Ephemeroptera_|_Plecoptera_Insecta_|_Trichoptera_|_Chironomidae", names(df.2_notnull))]

# Trim colnames to just the species
names <- gsub("^(\\w+_\\w+)_\\w+_\\w+_*\\w*", "\\1", colnames(df.2_notnull))
names(df.2_notnull)<-names

# move rownames to first col
setDT(df.2_notnull, keep.rownames = TRUE)[]

# put primer into own column
df.2_notnull$primer <- str_sub(df.2_notnull$rn, -1,-1)

# add a column of totals
species.df<-adorn_totals(df.2_notnull, where="col")

# put into summary table for ggplot
species.df2<-data.frame(species.df$Total, species.df$primer)

#rename cols
names(species.df2) <- c("pooled","primer")

# add column for rank
species.df2$rank<-"species"

# add column for indicator
species.df2$ind<-"fresh"

#######################################################
# Create dataframe for Genus rarefaction
######################################################

# Split up SampleName for Arthropoda matrix
H2<-data.frame(H, do.call(rbind, str_split(H$SampleName,"_")))

# Add new column names
names(H2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
H2_A<-H2[grep("^A",H2$A_GlobalESV),]
H2_A$primer<-"A"
H2_B<-H2[grep("^B",H2$A_GlobalESV),]
H2_B$primer<-"B"
H2_C<-H2[grep("^C",H2$A_GlobalESV),]
H2_C$primer<-"C"
H2_D<-H2[grep("^D",H2$A_GlobalESV),]
H2_D$primer<-"D"
H2_E<-H2[grep("^E",H2$A_GlobalESV),]
H2_E$primer<-"E"
H2_F<-H2[grep("^F",H2$A_GlobalESV),]
H2_F$primer<-"F"

#combine them all into a single dataframe
H3<-rbind(H2_A, H2_B, H2_C, H2_D, H2_E, H2_F)

# pivot to make Genus matrix only (pool data across sites and reps)
genus<-dcast(H3, kit+primer ~ Genus+Order+Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+primer into single column
genus$sample<-paste(genus$kit, genus$primer, sep="_")

# delete duplicate columns
genus$kit<-NULL
genus$primer<-NULL

#move SampleName to rownames then delete
rownames(genus)<-genus$sample
genus$sample<-NULL

#remove columns with only zeros
genus_notnull<-genus[,colSums(genus) !=0]

#remove rows with only zeros
genus_notnull2<-genus_notnull[rowSums(genus_notnull) !=0,]

#calculate 15th percentile for rrarefy function
genus_15percentile<-quantile(rowSums(genus_notnull2), prob=0.15)

# Normalize dataset down to 15th percentile
set.seed(1234)
df<-rrarefy(genus_notnull2, sample=genus_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset  the markers, work with soil kit only
df.2<-df.2[grepl("^Soil", rownames(df.2)),]

# Remove zero sum columns
df.2_notnull<-df.2[,colSums(df.2) !=0]

# Retain freshwater groups
df.2_notnull<-df.2_notnull[,grepl("_Ephemeroptera_|_Plecoptera_Insecta_|_Trichoptera_|_Chironomidae", names(df.2_notnull))]

# Trim colnames to just the genus
names <- gsub("^(\\w+)_\\w+_\\w+_*\\w*", "\\1", colnames(df.2_notnull))
names(df.2_notnull)<-names

# move rownames to first col
setDT(df.2_notnull, keep.rownames = TRUE)[]

# put primer into own column
df.2_notnull$primer <- str_sub(df.2_notnull$rn, -1,-1)

# add a column of totals
genus.df<-adorn_totals(df.2_notnull, where="col")

# put into summary table for ggplot
genus.df2<-data.frame(genus.df$Total, genus.df$primer)

#rename cols
names(genus.df2) <- c("pooled","primer")

# add column for rank
genus.df2$rank<-"genus"

# add column for indicator
genus.df2$ind<-"fresh"

#######################################################
# Create dataframe for Family rarefaction
######################################################

# pivot to make Family matrix only (pool data across sites and reps)
family<-dcast(E3, kit+primer ~ Family+Order+Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+primer into single column
family$sample<-paste(family$kit, family$primer, sep="_")

# delete duplicate columns
family$kit<-NULL
family$primer<-NULL

#move SampleName to rownames then delete
rownames(family)<-family$sample
family$sample<-NULL

#remove columns with only zeros
family_notnull<-family[,colSums(family) !=0]

#remove rows with only zeros & edit rownames
family_notnull2<-family_notnull[rowSums(family_notnull) !=0,]

#calculate 15th percentile for rrarefy function
family_15percentile<-quantile(rowSums(family_notnull2), prob=0.15)

# Normalize dataset down to 15th percentile
set.seed(1234)
df<-rrarefy(family_notnull2, sample=family_15percentile)

# Convert to presence-absence matrix
df[df>0]<-1

# Convert from type 'double' to 'list'
df.2 <-data.frame(df)

# Subset  the markers, work with soil only
df.2<-df.2[grepl("^Soil", rownames(df.2)),]

# Remove zero sum columns
df.2_notnull<-df.2[,colSums(df.2) !=0]

# Retain freshwater groups
df.2_notnull<-df.2_notnull[,grepl("_Ephemeroptera_|_Plecoptera_Insecta_|_Trichoptera_|_Chironomidae", names(df.2_notnull))]

# Trim colnames to just the family
names <- gsub("^(\\w+)_\\w+_\\w+_*\\w*", "\\1", colnames(df.2_notnull))
names(df.2_notnull)<-names

# move rownames to first col
setDT(df.2_notnull, keep.rownames = TRUE)[]

# put primer into own column
df.2_notnull$primer <- str_sub(df.2_notnull$rn, -1,-1)

# add a column of totals
family.df<-adorn_totals(df.2_notnull, where="col")

# put into summary table for ggplot
family.df2<-data.frame(family.df$Total, family.df$primer)

#rename cols
names(family.df2) <- c("pooled","primer")

# add column for rank
family.df2$rank<-"family"

# add column for indicator
family.df2$ind<-"fresh"

# Put all the summry results in ONE table
summary2<-rbind(esv.df2, species.df2, genus.df2, family.df2)

# Combine the broadscale and freshwater summaries
summary3<-rbind(summary, summary2)

# Create factors
summary3$rank<-factor(summary3$rank, levels=c("esv","species","genus","family"),
                      labels=c("ESV","Species","Genus","Family"))
summary3$primer<-factor(summary3$primer, levels=c("A","B","C","D","E","F"),
                        labels=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"))
summary3$ind<-factor(summary3$ind, levels=c("broad","fresh"), labels=c("Arthropoda Site Indicators","EPTC"))

# Create summary bar chart with facet wrap (rank)
custom<-brewer.pal(n = 7, name = 'Set1')
custom2<-custom[-6]

# Create a df to plot a maximum value point that is 'hidden'
temp <- data.frame(
                   ind=c(rep("EPTC",24),rep("Arthropoda Site Indicators",24)),
                   rank = c(rep("ESV",6),rep("Species",6),rep("Genus",6),rep("Family",6)),
                   primer=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"),
                   value = c(rep(520,6),rep(70,6),rep(70,6),rep(70,6),
                             rep(70,6),rep(20,6),rep(20,6),rep(20,6)))

# Create a facet wrap plot
ind<-ggplot(data=summary3, aes(x=primer, y=pooled, fill=primer)) +
  geom_bar(data=summary3, stat="identity") +
  geom_point(data = temp, aes(x = primer, y = value), colour = "white") +
  labs(x="Primer", y="Indicators") +
  facet_wrap(~ind+rank, scales="free", ncol=4) +
  scale_fill_manual(values=custom2) +
  scale_y_continuous(labels=comma) +
  theme(legend.title=element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title=element_text(size=10), 
        axis.text.x = element_text(angle=90))+
  guides(fill = guide_legend(nrow = 1))

ggsave("F3_site_waterquality_indicators.pdf", ind)
