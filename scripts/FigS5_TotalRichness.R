# Teresita M. Porter, May 30, 2019

library(stringr)
library(scales)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(ggplot2)

#####################################################################
# Read in sample x taxonomy table use this for all taxa
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Create matrix for just the arthropoda (use for ESVs, orders, classes)
B<-A[A$Phylum=="Arthropoda",]

# Create matrix for Arthropoda species with sBP>=0.70 (95% correct)
C<-B[B$sBP>=0.70,]

# Create matrix for Arthropoda genera with gBP>=0.30
D<-B[B$gBP>=0.30,]

# Create matrix for Arthropoda families with fBP>=0.20
E<-B[B$fBP>=0.20,]

#######################################################
# Create dataframe for rarefaction curves
## pool data across reps

# Split up SampleName for each matrix
A2<-data.frame(A, do.call(rbind, str_split(A$SampleName,"_")))
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
C2<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")))
D2<-data.frame(D, do.call(rbind, str_split(D$SampleName,"_")))
E2<-data.frame(E, do.call(rbind, str_split(E$SampleName,"_")))

# Add new column names
names(A2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(B2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(C2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(D2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(E2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

# create df for each primer and add primer name to a new primer column
A2_A<-A2[grep("^A",A2$A_GlobalESV),]
A2_A$primer<-"A"
A2_B<-A2[grep("^B",A2$A_GlobalESV),]
A2_B$primer<-"B"
A2_C<-A2[grep("^C",A2$A_GlobalESV),]
A2_C$primer<-"C"
A2_D<-A2[grep("^D",A2$A_GlobalESV),]
A2_D$primer<-"D"
A2_E<-A2[grep("^E",A2$A_GlobalESV),]
A2_E$primer<-"E"
A2_F<-A2[grep("^F",A2$A_GlobalESV),]
A2_F$primer<-"F"
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
A3<-rbind(A2_A, A2_B, A2_C, A2_D, A2_E, A2_F)
B3<-rbind(B2_A, B2_B, B2_C, B2_D, B2_E, B2_F)
C3<-rbind(C2_A, C2_B, C2_C, C2_D, C2_E, C2_F)
D3<-rbind(D2_A, D2_B, D2_C, D2_D, D2_E, D2_F)
E3<-rbind(E2_A, E2_B, E2_C, E2_D, E2_E, E2_F)

# pivot to make esv matrix only
all_esv<-dcast(A3, kit+primer+site ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
esv<-dcast(B3, kit+primer+site ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
species<-dcast(C3, kit+primer+site ~ Species, value.var = "ESVsize", fun.aggregate = sum)
genus<-dcast(D3, kit+primer+site ~ Genus, value.var = "ESVsize", fun.aggregate = sum)
family<-dcast(E3, kit+primer+site ~ Family, value.var = "ESVsize", fun.aggregate = sum)
order<-dcast(A3, kit+primer+site ~ Order, value.var = "ESVsize", fun.aggregate = sum)
class<-dcast(A3, kit+primer+site ~ Class, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
all_esv$sample<-paste(all_esv$kit, all_esv$primer, all_esv$site, sep="_")
esv$sample<-paste(esv$kit, esv$primer, esv$site, sep="_")
species$sample<-paste(species$kit, species$primer, species$site, sep="_")
genus$sample<-paste(genus$kit, genus$primer, genus$site, sep="_")
family$sample<-paste(family$kit, family$primer, family$site, sep="_")
order$sample<-paste(order$kit, order$primer, order$site, sep="_")
class$sample<-paste(class$kit, class$primer, class$site, sep="_")

# delete duplicate columns
all_esv$kit<-NULL
all_esv$primer<-NULL
all_esv$site<-NULL
esv$kit<-NULL
esv$primer<-NULL
esv$site<-NULL
species$kit<-NULL
species$primer<-NULL
species$site<-NULL
genus$kit<-NULL
genus$primer<-NULL
genus$site<-NULL
family$kit<-NULL
family$primer<-NULL
family$site<-NULL
order$kit<-NULL
order$primer<-NULL
order$site<-NULL
class$kit<-NULL
class$primer<-NULL
class$site<-NULL

#move SampleName to rownames then delete
rownames(all_esv)<-all_esv$sample
all_esv$sample<-NULL
rownames(esv)<-esv$sample
esv$sample<-NULL
rownames(species)<-species$sample
species$sample<-NULL
rownames(genus)<-genus$sample
genus$sample<-NULL
rownames(family)<-family$sample
family$sample<-NULL
rownames(order)<-order$sample
order$sample<-NULL
rownames(class)<-class$sample
class$sample<-NULL

#remove columns with only zeross
all_esv_notnull<-all_esv[,colSums(all_esv) !=0]
esv_notnull<-esv[,colSums(esv) !=0]
species_notnull<-species[,colSums(species) !=0]
genus_notnull<-genus[,colSums(genus) !=0]
family_notnull<-family[,colSums(family) !=0]
order_notnull<-order[,colSums(order) !=0]
class_notnull<-class[,colSums(class) !=0]

#remove rows with only zeros & edit rownames
all_esv_notnull2<-all_esv_notnull[rowSums(all_esv_notnull) !=0,]
esv_notnull2<-esv_notnull[rowSums(esv_notnull) !=0,]
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]
genus_notnull2<-genus_notnull[rowSums(genus_notnull) !=0,]
family_notnull2<-family_notnull[rowSums(family_notnull) !=0,]
order_notnull2<-order_notnull[rowSums(order_notnull) !=0,]
class_notnull2<-class_notnull[rowSums(class_notnull) !=0,]

#calculate 15th percentile for rrarefy function
all_esv_15percentile<-quantile(rowSums(all_esv_notnull2), prob=0.15)
esv_15percentile<-quantile(rowSums(esv_notnull2), prob=0.15)
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)
genus_15percentile<-quantile(rowSums(genus_notnull2), prob=0.15)
family_15percentile<-quantile(rowSums(family_notnull2), prob=0.15)
order_15percentile<-quantile(rowSums(order_notnull2), prob=0.15)
class_15percentile<-quantile(rowSums(class_notnull2), prob=0.15)

###################################################################
##### Rarefy dataset
###################################################################

set.seed(1234)

#Rarefy original ESV matrix down to 15th percentile library size to normalize read depth across samples

all_esv_df<-rrarefy(all_esv_notnull2, sample=all_esv_15percentile)
esv_df<-rrarefy(esv_notnull2, sample=esv_15percentile)
species_df<-rrarefy(species_notnull2, sample=species_15percentile)
genus_df<-rrarefy(genus_notnull2, sample=genus_15percentile)
family_df<-rrarefy(family_notnull2, sample=family_15percentile)
order_df<-rrarefy(order_notnull2, sample=order_15percentile)
class_df<-rrarefy(class_notnull2, sample=class_15percentile)

###

# Process all_esvs get uniques SOIL
# Count all unique ESVs then add to dataframe for ggplot cols uniques, rank
# subset
temp<-all_esv_df[grep("^SoilKit",rownames(all_esv_df)),]
#drop columns with all zeros
temp<-temp[,colSums(temp) !=0]
uniques<-data.frame(length(unique(colnames(temp))), "ESV","all","all","soil")
names(uniques)<-c("uniques","rank","primer","taxa","kit")
## for A
temp<-all_esv_df[grep("^SoilKit_A_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","A","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-all_esv_df[grep("^SoilKit_B_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","B","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-all_esv_df[grep("^SoilKit_C_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","C","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-all_esv_df[grep("^SoilKit_D_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","D","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-all_esv_df[grep("^SoilKit_E_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","E","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-all_esv_df[grep("^SoilKit_F_",rownames(all_esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","F","all","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda esvs get uniques
# Count all unique ESVs then add to dataframe for ggplot cols uniques, rank
temp<-esv_df[grep("^SoilKit_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-esv_df[grep("^SoilKit_A_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-esv_df[grep("^SoilKit_B_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-esv_df[grep("^SoilKit_C_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-esv_df[grep("^SoilKit_D_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-esv_df[grep("^SoilKit_E_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-esv_df[grep("^SoilKit_F_",rownames(esv_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "ESV","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda species get uniques
# Count all unique species then add to dataframe for ggplot cols uniques, rank
temp<-species_df[grep("^SoilKit_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-species_df[grep("^SoilKit_A_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-species_df[grep("^SoilKit_B_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-species_df[grep("^SoilKit_C_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-species_df[grep("^SoilKit_D_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-species_df[grep("^SoilKit_E_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-species_df[grep("^SoilKit_F_",rownames(species_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "species","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda genera get uniques
# Count all unique genera then add to dataframe for ggplot cols uniques, rank
temp<-genus_df[grep("^SoilKit_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-genus_df[grep("^SoilKit_A_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-genus_df[grep("^SoilKit_B_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-genus_df[grep("^SoilKit_C_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-genus_df[grep("^SoilKit_D_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-genus_df[grep("^SoilKit_E_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-genus_df[grep("^SoilKit_F_",rownames(genus_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "genus","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda families get uniques
# Count all unique families then add to dataframe for ggplot cols uniques, rank
temp<-family_df[grep("^SoilKit_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-family_df[grep("^SoilKit_A_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-family_df[grep("^SoilKit_B_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-family_df[grep("^SoilKit_C_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-family_df[grep("^SoilKit_D_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-family_df[grep("^SoilKit_E_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-family_df[grep("^SoilKit_F_",rownames(family_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "family","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda orders get uniques
# Count all unique orders then add to dataframe for ggplot cols uniques, rank
temp<-order_df[grep("^SoilKit_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-order_df[grep("^SoilKit_A_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-order_df[grep("^SoilKit_B_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-order_df[grep("^SoilKit_C_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-order_df[grep("^SoilKit_D_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-order_df[grep("^SoilKit_E_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-order_df[grep("^SoilKit_F_",rownames(order_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "order","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# Process arthropoda classes get uniques
# Count all unique classes then add to dataframe for ggplot cols uniques, rank
temp<-class_df[grep("^SoilKit_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","all","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for A
temp<-class_df[grep("^SoilKit_A_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","A","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for B
temp<-class_df[grep("^SoilKit_B_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","B","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for C
temp<-class_df[grep("^SoilKit_C_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","C","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for D
temp<-class_df[grep("^SoilKit_D_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","D","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for E
temp<-class_df[grep("^SoilKit_E_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","E","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)
## for F
temp<-class_df[grep("^SoilKit_F_",rownames(class_df)),]
temp<-temp[,colSums(temp) !=0]
uniques.tmp<-data.frame(length(unique(colnames(temp))), "class","F","arthropoda","soil")
names(uniques.tmp)<-c("uniques","rank","primer","taxa","kit")
uniques<-rbind(uniques,uniques.tmp)

# merge taxa with rank into single column
uniques$kit_taxa_rank <- paste(uniques$kit, uniques$taxa, uniques$rank, sep="_")

# Create factors
uniques$kit_taxa_rank <- factor(uniques$kit_taxa_rank, 
                            levels=c("soil_all_ESV","soil_arthropoda_ESV","soil_arthropoda_species","soil_arthropoda_genus","soil_arthropoda_family","soil_arthropoda_order","soil_arthropoda_class","tissue_all_ESV","tissue_arthropoda_ESV","tissue_arthropoda_species","tissue_arthropoda_genus","tissue_arthropoda_family","tissue_arthropoda_order","tissue_arthropoda_class"),
                            labels=c("(S) All ESVs","(S) ESVs","(S) Species","(S) Genera","(S) Families","(S) Orders","(S) Classes","(T) All ESVs","(T) ESVs","(T) Species","(T) Genera","(T) Families","(T) Orders","(T) Classes"))

uniques$primer <- factor(uniques$primer,
                         levels=c("all","A","B","C","D","E","F"),
                         labels=c("All","BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"))

uniques$kit<-factor(uniques$kit,
                    levels=c("soil","tissue"),
                    labels=c("Soil","Tissue"))
uniques$rank<-factor(uniques$rank,
                     levels=c("ESV","species","genus","family","order","class"),
                     labels=c("ESV","Species","Genus","Family","Order","Class"))

# Only plot data for arthropoda (represents ~ 25% of all unique ESVs)
arth_uniques<-uniques[uniques$taxa!="all",]

###########################
# Describe arthropoda results
###########################

custom<-brewer.pal(n = 7, name = 'Set1')
custom2<-custom[-6]
custom3<-c('#999999',custom2)

p1<-ggplot(data=arth_uniques, aes(x=primer,y=uniques, fill=primer)) +
  geom_bar(stat="identity") +
  labs(x="Primers", y="Richness") +
  scale_y_continuous(labels=comma) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.title=element_blank()) +
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~rank, scales="free", ncol=6) +
  scale_fill_manual(values=custom3)

ggsave("SF5_TotalRichness.pdf",p1)

# Visual tests for normality (also tried for each kit separately and pooled)
ggdensity(arth_uniques$uniques[arth_uniques$kit=="Soil"], 
          main = "density uniques dist",
          xlab = "uniques")

ggqqplot(arth_uniques$uniques[arth_uniques$kit=="Soil"])

# Shapiro-Wilk test of normality (often positive for small sample sizes)
shapiro.test(arth_uniques$uniques)
# W = 0.38107, p-value = 1.751e-14 show significant difference from normal distribution too

# Kruskal-Wallis one way analysis of variance
kruskal.test(primer ~ rank, data = arth_uniques[(arth_uniques$kit=="Soil") & (arth_uniques$rank=="ESV"),]) 
