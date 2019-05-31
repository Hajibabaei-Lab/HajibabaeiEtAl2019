# Teresita M. Porter, May 30, 2019

library(stringr)
library(scales)
library(reshape2)
library(RColorBrewer)
library(vegan)
library(ggplot2)
library(gridExtra)

########################################################
# NEW FUNCTION TO CHOOSE TWO MARKERS AT A TIME
# df should be the rank_df
# marker1 and marker2 from combn function
choose_two<-function(df, marker1, marker2){
  
  p1<-paste("^", marker1, "_", sep="")
  p2<-paste("^", marker2, "_", sep="")
  
  target<-rownames(df)
  
  # grep the rownames
  df2<-df[( (grepl(p1, target)) |
              (grepl(p2, target)) ),]
  
  newList<-list("matrix"=df2, "pattern1"=p1, "pattern2"=p2, "target"=target)
  return(newList)
  
}
#####################################################################
# NEW FUNCTION TO CHOOSE THREE MARKERS AT A TIME
# df should be the rank_df
# marker1, marker2, and marker3 from combn function
choose_three<-function(df, marker1, marker2, marker3){
  
  p1<-paste("^", marker1, "_", sep="")
  p2<-paste("^", marker2, "_", sep="")
  p3<-paste("^", marker3, "_", sep="")
  
  target<-rownames(df)
  
  # grep the rownames
  df2<-df[( (grepl(p1, target)) |
              (grepl(p2, target)) |
              (grepl(p3, target)) ),]
  
  newList<-list("matrix"=df2, "pattern1"=p1, "pattern2"=p2, "pattern3"=p3, "target"=target)
  return(newList)
  
}

#####################################################################
# NEW FUNCTION TO CHOOSE FOUR MARKERS AT A TIME
# df should be the rank_df
# marker1, marker2, marker3, and marker4 from combn function
choose_four<-function(df, marker1, marker2, marker3, marker4){
  
  p1<-paste("^", marker1, "_", sep="")
  p2<-paste("^", marker2, "_", sep="")
  p3<-paste("^", marker3, "_", sep="")
  p4<-paste("^", marker4, "_", sep="")
  
  target<-rownames(df)
  
  # grep the rownames
  df2<-df[( (grepl(p1, target)) |
              (grepl(p2, target)) |
              (grepl(p3, target)) |
              (grepl(p4, target)) ),]
  
  newList<-list("matrix"=df2, "pattern1"=p1, "pattern2"=p2, "pattern3"=p3, "pattern4"=p4, "target"=target)
  return(newList)
  
}
#####################################################################
# NEW FUNCTION TO CHOOSE FIVE MARKERS AT A TIME
# df should be the rank_df
# marker1, marker2, marker3, marker4, and marker5 from combn function
choose_five<-function(df, marker1, marker2, marker3, marker4, marker5){
  
  p1<-paste("^", marker1, "_", sep="")
  p2<-paste("^", marker2, "_", sep="")
  p3<-paste("^", marker3, "_", sep="")
  p4<-paste("^", marker4, "_", sep="")
  p5<-paste("^", marker5, "_", sep="")
  
  target<-rownames(df)
  
  # grep the rownames
  df2<-df[( (grepl(p1, target)) |
              (grepl(p2, target)) |
              (grepl(p3, target)) |
              (grepl(p4, target)) |
              (grepl(p5, target)) ),]
  
  newList<-list("matrix"=df2, "pattern1"=p1, "pattern2"=p2, "pattern3"=p3, "pattern4"=p4, "pattern5"=p5, "target"=target)
  return(newList)
  
}
#####################################################################
# NEW FUNCTION TO CHOOSE SIX MARKERS AT A TIME
# df should be the rank_df
# marker1, marker2, marker3, marker4, marker5 and marker 6 from combn function
choose_six<-function(df, marker1, marker2, marker3, marker4, marker5, marker6){
  
  p1<-paste("^", marker1, "_", sep="")
  p2<-paste("^", marker2, "_", sep="")
  p3<-paste("^", marker3, "_", sep="")
  p4<-paste("^", marker4, "_", sep="")
  p5<-paste("^", marker5, "_", sep="")
  p6<-paste("^", marker6, "_", sep="")
  
  target<-rownames(df)
  
  # grep the rownames
  df2<-df[( (grepl(p1, target)) |
              (grepl(p2, target)) |
              (grepl(p3, target)) |
              (grepl(p4, target)) |
              (grepl(p5, target)) |
              (grepl(p6, target)) ),]
  
  newList<-list("matrix"=df2, "pattern1"=p1, "pattern2"=p2, "pattern3"=p3, "pattern4"=p4, "pattern5"=p5, "pattern6"=p6, "target"=target)
  return(newList)
  
}
#####################################################################

# Read in sample x taxonomy table use this for all taxa
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Create matrix for just the arthropoda 
B<-A[A$Phylum=="Arthropoda",]

# Work with soil kit only
B<-B[grepl("SoilKit", B$SampleName),] #(also use for ESVs, orders)

# Create matrix for Arthropoda species with sBP>=0.70 (95% correct)
C<-B[B$sBP>=0.70,]

# Create matrix for Arthropoda genera with gBP>=0.30 (99% correct)
D<-B[B$gBP>=0.30,]

# Create matrix for Arthropoda families with fBP>=0.20 (99% correct)
E<-B[B$fBP>=0.20,]

#######################################################
# Create dataframe for rarefaction curves
## pool data across reps

# Split up SampleName for each matrix
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))
C2<-data.frame(C, do.call(rbind, str_split(C$SampleName,"_")))
D2<-data.frame(D, do.call(rbind, str_split(D$SampleName,"_")))
E2<-data.frame(E, do.call(rbind, str_split(E$SampleName,"_")))

# Add new column names
names(B2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(C2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(D2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")
names(E2)[32:39]<-c("project","kit","site","rep","marker1","marker2","run","lane")

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

# Combine them all into a single dataframe
B3<-rbind(B2_A, B2_B, B2_C, B2_D, B2_E, B2_F)
C3<-rbind(C2_A, C2_B, C2_C, C2_D, C2_E, C2_F)
D3<-rbind(D2_A, D2_B, D2_C, D2_D, D2_E, D2_F)
E3<-rbind(E2_A, E2_B, E2_C, E2_D, E2_E, E2_F)

# Pivot to make rank matrices only
#all_ESV<-dcast(A3, primer+site ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
ESV<-dcast(B3, primer+site ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
species<-dcast(C3, primer+site ~ Species, value.var = "ESVsize", fun.aggregate = sum)
genus<-dcast(D3, primer+site ~ Genus, value.var = "ESVsize", fun.aggregate = sum)
family<-dcast(E3, primer+site ~ Family, value.var = "ESVsize", fun.aggregate = sum)
order<-dcast(B3, primer+site ~ Order, value.var = "ESVsize", fun.aggregate = sum)
class<-dcast(B3, primer+site ~ Class, value.var = "ESVsize", fun.aggregate = sum)

# Merge kit+site+rep+primer into single column
ESV$sample<-paste(ESV$primer, ESV$site, sep="_")
species$sample<-paste(species$primer, species$site, sep="_")
genus$sample<-paste(genus$primer, genus$site, sep="_")
family$sample<-paste(family$primer, family$site, sep="_")
order$sample<-paste(order$primer, order$site, sep="_")
class$sample<-paste(class$primer, class$site, sep="_")

# Delete duplicate columns
ESV$primer<-NULL
ESV$site<-NULL
species$primer<-NULL
species$site<-NULL
genus$primer<-NULL
genus$site<-NULL
family$primer<-NULL
family$site<-NULL
order$primer<-NULL
order$site<-NULL
class$primer<-NULL
class$site<-NULL

# Move SampleName to rownames then delete
rownames(ESV)<-ESV$sample
ESV$sample<-NULL
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
ESV_notnull<-ESV[,colSums(ESV) !=0]
species_notnull<-species[,colSums(species) !=0]
genus_notnull<-genus[,colSums(genus) !=0]
family_notnull<-family[,colSums(family) !=0]
order_notnull<-order[,colSums(order) !=0]
class_notnull<-class[,colSums(class) !=0]

#remove rows with only zeros & edit rownames
ESV_notnull2<-ESV_notnull[rowSums(ESV_notnull) !=0,]
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]
genus_notnull2<-genus_notnull[rowSums(genus_notnull) !=0,]
family_notnull2<-family_notnull[rowSums(family_notnull) !=0,]
order_notnull2<-order_notnull[rowSums(order_notnull) !=0,]
class_notnull2<-class_notnull[rowSums(class_notnull) !=0,]

# Calculate 15th percentile for rrarefy function
ESV_15percentile<-quantile(rowSums(ESV_notnull2), prob=0.15)
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)
genus_15percentile<-quantile(rowSums(genus_notnull2), prob=0.15)
family_15percentile<-quantile(rowSums(family_notnull2), prob=0.15)
order_15percentile<-quantile(rowSums(order_notnull2), prob=0.15)
class_15percentile<-quantile(rowSums(class_notnull2), prob=0.15)

###################################################################
##### Rarefy dataset
###################################################################

set.seed(1234)
#Rarefy original species matrix down to 15th percentile library size to normalize read depth across samples

esv_df<-rrarefy(ESV_notnull2, sample=ESV_15percentile)
species_df<-rrarefy(species_notnull2, sample=species_15percentile)
genus_df<-rrarefy(genus_notnull2, sample=genus_15percentile)
family_df<-rrarefy(family_notnull2, sample=family_15percentile)
order_df<-rrarefy(order_notnull2, sample=order_15percentile)
class_df<-rrarefy(class_notnull2, sample=class_15percentile)

# Create list of markers
markers=c("A","B","C","D","E","F")

##############################################
# Create tables for ESV rank
##############################################

# initialize single primer results dataframe
single_results<-data.frame(marker=character(),
                           uniques=integer())

# get number of unique taxa per marker
for (i in markers) {
  p1<-paste("^",i,"_",sep="")
  target<-rownames(esv_df)
  x<-esv_df[grep(p1,target),]
  y<-x[,colSums(x) !=0]
  single_results.tmp<-data.frame(marker=i, uniques=length(unique(colnames(y))))
  single_results<-rbind(single_results,single_results.tmp)
}

# sorted by descending number of uniques
single_sorted_esv<-single_results[order(-single_results$uniques),]

# add column to indicate number of markers combined
single_sorted_esv$combo<-"1"

# choose markers in pairwise fashion
pairwise<-data.frame(combn(markers, 2))

# initialize results dataframe
double_results<-data.frame(marker1=character(),
                    marker2=character(),
                    uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(pairwise)) {
  marker1<-pairwise[[i]][1]
  marker2<-pairwise[[i]][2]
  test<-choose_two(esv_df, marker1, marker2)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  double_results.tmp<-data.frame(marker1=marker1, marker2=marker2, uniques=length(unique(colnames(y))))
  double_results<-rbind(double_results,double_results.tmp)
}

# sorted by descending number of uniques
double_sorted_esv<-double_results[order(-double_results$uniques),]

#reformat then merge with single_sorted_esv
double_sorted_esv$marker<-paste(double_sorted_esv$marker1, double_sorted_esv$marker2)
double_sorted_esv<-double_sorted_esv[,-c(1:2)]
double_sorted_esv$combo<-"2"
esv<-rbind(single_sorted_esv, double_sorted_esv)

# Repeat for 3 markers
# Chooses markers in three-way fashion
threeway<-data.frame(combn(markers, 3))

# initialize results dataframe
triple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(threeway)) {
  marker1<-threeway[[i]][1]
  marker2<-threeway[[i]][2]
  marker3<-threeway[[i]][3]
  test<-choose_three(esv_df, marker1, marker2, marker3)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  triple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, uniques=length(unique(colnames(y))))
  triple_results<-rbind(triple_results,triple_results.tmp)
}

# sorted by descending number of uniques
triple_sorted_esv<-triple_results[order(-triple_results$uniques),]

#reformat then merge with esv
triple_sorted_esv$marker<-paste(triple_sorted_esv$marker1, triple_sorted_esv$marker2, triple_sorted_esv$marker3)
triple_sorted_esv<-triple_sorted_esv[,-c(1:3)]
triple_sorted_esv$combo<-"3"
esv<-rbind(esv, triple_sorted_esv)

# Repeat for 4 markers
# Chooses markers in four-way fashion
fourway<-data.frame(combn(markers, 4))

# initialize results dataframe
quadruple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           marker4=character(),
                           uniques=integer())



## for each column in pairwise, run combo_uniques
for (i in names(fourway)) {
  marker1<-fourway[[i]][1]
  marker2<-fourway[[i]][2]
  marker3<-fourway[[i]][3]
  marker4<-fourway[[i]][4]
  test<-choose_four(esv_df, marker1, marker2, marker3, marker4)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quadruple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, uniques=length(unique(colnames(y))))
  quadruple_results<-rbind(quadruple_results,quadruple_results.tmp)
}

# sorted by descending number of uniques
quadruple_sorted_esv<-quadruple_results[order(-quadruple_results$uniques),]

#reformat then merge with esv
quadruple_sorted_esv$marker<-paste(quadruple_sorted_esv$marker1, quadruple_sorted_esv$marker2, quadruple_sorted_esv$marker3, quadruple_sorted_esv$marker4)
quadruple_sorted_esv<-quadruple_sorted_esv[,-c(1:4)]
quadruple_sorted_esv$combo<-"4"
esv<-rbind(esv, quadruple_sorted_esv)

# Repeat for 5 markers
# Chooses markers in four-way fashion
fiveway<-data.frame(combn(markers, 5))

# initialize results dataframe
quintuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fiveway)) {
  marker1<-fiveway[[i]][1]
  marker2<-fiveway[[i]][2]
  marker3<-fiveway[[i]][3]
  marker4<-fiveway[[i]][4]
  marker5<-fiveway[[i]][5]
  
  test<-choose_five(esv_df, marker1, marker2, marker3, marker4, marker5)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quintuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, uniques=length(unique(colnames(y))))
  quintuple_results<-rbind(quintuple_results,quintuple_results.tmp)
}

# sorted by descending number of uniques
quintuple_sorted_esv<-quintuple_results[order(-quintuple_results$uniques),]

#reformat then merge with esv
quintuple_sorted_esv$marker<-paste(quintuple_sorted_esv$marker1, quintuple_sorted_esv$marker2, quintuple_sorted_esv$marker3, quintuple_sorted_esv$marker4, quintuple_sorted_esv$marker5)
quintuple_sorted_esv<-quintuple_sorted_esv[,-c(1:5)]
quintuple_sorted_esv$combo<-"5"
esv<-rbind(esv, quintuple_sorted_esv)

# Repeat for 6 markers
# Chooses markers in six-way fashion
sixway<-data.frame(combn(markers, 6))

# initialize results dataframe
sextuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              marker6=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(sixway)) {
  marker1<-sixway[[i]][1]
  marker2<-sixway[[i]][2]
  marker3<-sixway[[i]][3]
  marker4<-sixway[[i]][4]
  marker5<-sixway[[i]][5]
  marker6<-sixway[[i]][6]
  
  test<-choose_six(esv_df, marker1, marker2, marker3, marker4, marker5, marker6)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  sextuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, marker6=marker6, uniques=length(unique(colnames(y))))
  sextuple_results<-rbind(sextuple_results, sextuple_results.tmp)
}

# sorted by descending number of uniques
sextuple_sorted_esv<-sextuple_results[order(-sextuple_results$uniques),]

#reformat then merge with esv
sextuple_sorted_esv$marker<-paste(sextuple_sorted_esv$marker1, sextuple_sorted_esv$marker2, sextuple_sorted_esv$marker3, sextuple_sorted_esv$marker4, sextuple_sorted_esv$marker5, sextuple_sorted_esv$marker6)
sextuple_sorted_esv<-sextuple_sorted_esv[,-c(1:6)]
sextuple_sorted_esv$combo<-"6"
esv<-rbind(esv, sextuple_sorted_esv)
#create copy for editing
esv2<-esv

# replace letter codes with primer names
esv2$marker<-gsub("F", "fwh1", esv2$marker)
esv2$marker<-gsub("B", "F230R", esv2$marker)
esv2$marker<-gsub("A","BR5",esv2$marker)
esv2$marker<-gsub("C", "ml-jg", esv2$marker)
esv2$marker<-gsub("D", "BF1R2", esv2$marker)
esv2$marker<-gsub("E", "BF2R2", esv2$marker)

# Plot an accumulation curve
f1<-ggplot(esv2, aes(x=combo, y=uniques, label = marker)) +
  ggtitle("1)") +
  geom_point() +
  annotate("text", label="C\n", x=1, y=750, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="B + C\n", x=2, y=1350, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C /\n B + C + D", x=3, y=2330, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D", x=4, y=3030, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E", x=5, y=3486, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E + F", x=6, y=3746, hjust=0.5, vjust=-0.5, size=2.25) +
  ylim(1,4500) +
  labs(x="No. Combined Markers", y="ESV richness") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

##############################################
# Create tables for species rank
##############################################

# initialize single primer results dataframe
single_results<-data.frame(marker=character(),
                           uniques=integer())

# get number of unique taxa per marker
for (i in markers) {
  p1<-paste("^",i,"_",sep="")
  target<-rownames(species_df)
  x<-species_df[grep(p1,target),]
  y<-x[,colSums(x) !=0]
  single_results.tmp<-data.frame(marker=i, uniques=length(unique(colnames(y))))
  single_results<-rbind(single_results,single_results.tmp)
}

# sorted by descending number of uniques
single_sorted_species<-single_results[order(-single_results$uniques),]

# add column to indicate number of markers combined
single_sorted_species$combo<-"1"

# initialize double results dataframe
double_results<-data.frame(marker1=character(),
                           marker2=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(pairwise)) {
  marker1<-pairwise[[i]][1]
  marker2<-pairwise[[i]][2]
  test<-choose_two(species_df, marker1, marker2)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  double_results.tmp<-data.frame(marker1=marker1, marker2=marker2, uniques=length(unique(colnames(y))))
  double_results<-rbind(double_results,double_results.tmp)
}

# sorted by descending number of uniques
double_sorted_species<-double_results[order(-double_results$uniques),]

#reformat then merge with single_sorted_species
double_sorted_species$marker<-paste(double_sorted_species$marker1, double_sorted_species$marker2)
double_sorted_species<-double_sorted_species[,-c(1:2)]
double_sorted_species$combo<-"2"
species<-rbind(single_sorted_species, double_sorted_species)

# initialize triple results dataframe
triple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(threeway)) {
  marker1<-threeway[[i]][1]
  marker2<-threeway[[i]][2]
  marker3<-threeway[[i]][3]
  test<-choose_three(species_df, marker1, marker2, marker3)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  triple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, uniques=length(unique(colnames(y))))
  triple_results<-rbind(triple_results,triple_results.tmp)
}

# sorted by descending number of uniques
triple_sorted_species<-triple_results[order(-triple_results$uniques),]

#reformat then merge with species
triple_sorted_species$marker<-paste(triple_sorted_species$marker1, triple_sorted_species$marker2, triple_sorted_species$marker3)
triple_sorted_species<-triple_sorted_species[,-c(1:3)]
triple_sorted_species$combo<-"3"
species<-rbind(species, triple_sorted_species)

# initialize quadruple results dataframe
quadruple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fourway)) {
  marker1<-fourway[[i]][1]
  marker2<-fourway[[i]][2]
  marker3<-fourway[[i]][3]
  marker4<-fourway[[i]][4]
  test<-choose_four(species_df, marker1, marker2, marker3, marker4)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quadruple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, uniques=length(unique(colnames(y))))
  quadruple_results<-rbind(quadruple_results,quadruple_results.tmp)
}

# sorted by descending number of uniques
quadruple_sorted_species<-quadruple_results[order(-quadruple_results$uniques),]

#reformat then merge with species
quadruple_sorted_species$marker<-paste(quadruple_sorted_species$marker1, quadruple_sorted_species$marker2, quadruple_sorted_species$marker3, quadruple_sorted_species$marker4)
quadruple_sorted_species<-quadruple_sorted_species[,-c(1:4)]
quadruple_sorted_species$combo<-"4"
species<-rbind(species, quadruple_sorted_species)

# initialize quintuple results dataframe
quintuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fiveway)) {
  marker1<-fiveway[[i]][1]
  marker2<-fiveway[[i]][2]
  marker3<-fiveway[[i]][3]
  marker4<-fiveway[[i]][4]
  marker5<-fiveway[[i]][5]
  
  test<-choose_five(species_df, marker1, marker2, marker3, marker4, marker5)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quintuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, uniques=length(unique(colnames(y))))
  quintuple_results<-rbind(quintuple_results,quintuple_results.tmp)
}

# sorted by descending number of uniques
quintuple_sorted_species<-quintuple_results[order(-quintuple_results$uniques),]

#reformat then merge with species
quintuple_sorted_species$marker<-paste(quintuple_sorted_species$marker1, quintuple_sorted_species$marker2, quintuple_sorted_species$marker3, quintuple_sorted_species$marker4, quintuple_sorted_species$marker5)
quintuple_sorted_species<-quintuple_sorted_species[,-c(1:5)]
quintuple_sorted_species$combo<-"5"
species<-rbind(species, quintuple_sorted_species)

# initialize sextuple results dataframe
sextuple_results<-data.frame(marker1=character(),
                             marker2=character(),
                             marker3=character(),
                             marker4=character(),
                             marker5=character(),
                             marker6=character(),
                             uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(sixway)) {
  marker1<-sixway[[i]][1]
  marker2<-sixway[[i]][2]
  marker3<-sixway[[i]][3]
  marker4<-sixway[[i]][4]
  marker5<-sixway[[i]][5]
  marker6<-sixway[[i]][6]
  
  test<-choose_six(species_df, marker1, marker2, marker3, marker4, marker5, marker6)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  sextuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, marker6=marker6, uniques=length(unique(colnames(y))))
  sextuple_results<-rbind(sextuple_results, sextuple_results.tmp)
}

# sorted by descending number of uniques
sextuple_sorted_species<-sextuple_results[order(-sextuple_results$uniques),]

#reformat then merge with species
sextuple_sorted_species$marker<-paste(sextuple_sorted_species$marker1, sextuple_sorted_species$marker2, sextuple_sorted_species$marker3, sextuple_sorted_species$marker4, sextuple_sorted_species$marker5, sextuple_sorted_species$marker6)
sextuple_sorted_species<-sextuple_sorted_species[,-c(1:6)]
sextuple_sorted_species$combo<-"6"
species<-rbind(species, sextuple_sorted_species)
#create copy for editing
species2<-species

# replace letter codes with primer names
species2$marker<-gsub("F", "fwh1", species2$marker)
species2$marker<-gsub("B", "F230R", species2$marker)
species2$marker<-gsub("A","BR5",species2$marker)
species2$marker<-gsub("C", "ml-jg", species2$marker)
species2$marker<-gsub("D", "BF1R2", species2$marker)
species2$marker<-gsub("E", "BF2R2", species2$marker)

# Plot an accumulation curve
f2<-ggplot(species2, aes(x=combo, y=uniques, label = marker)) +
  ggtitle("2)") +
  geom_point() +
  annotate("text", label="B\n", x=1, y=75, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B /\n B + C", x=2, y=96, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C", x=3, y=110, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D", x=4, y=111, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E", x=5, y=111, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E + F", x=6, y=111, hjust=0.5, vjust=-0.5, size=2.25) +
  labs(x="No. Combined Markers", y="Species richness") +
  ylim(1,150) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

##############################################
# Create tables for genus rank
##############################################

# initialize single primer results dataframe
single_results<-data.frame(marker=character(),
                           uniques=integer())

# get number of unique taxa per marker
for (i in markers) {
  p1<-paste("^",i,"_",sep="")
  target<-rownames(genus_df)
  x<-genus_df[grep(p1,target),]
  y<-x[,colSums(x) !=0]
  single_results.tmp<-data.frame(marker=i, uniques=length(unique(colnames(y))))
  single_results<-rbind(single_results,single_results.tmp)
}

# sorted by descending number of uniques
single_sorted_genus<-single_results[order(-single_results$uniques),]

# add column to indicate number of markers combined
single_sorted_genus$combo<-"1"

# initialize double results dataframe
double_results<-data.frame(marker1=character(),
                           marker2=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(pairwise)) {
  marker1<-pairwise[[i]][1]
  marker2<-pairwise[[i]][2]
  test<-choose_two(genus_df, marker1, marker2)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  double_results.tmp<-data.frame(marker1=marker1, marker2=marker2, uniques=length(unique(colnames(y))))
  double_results<-rbind(double_results,double_results.tmp)
}

# sorted by descending number of uniques
double_sorted_genus<-double_results[order(-double_results$uniques),]

#reformat then merge with single_sorted_genus
double_sorted_genus$marker<-paste(double_sorted_genus$marker1, double_sorted_genus$marker2)
double_sorted_genus<-double_sorted_genus[,-c(1:2)]
double_sorted_genus$combo<-"2"
genus<-rbind(single_sorted_genus, double_sorted_genus)

# initialize triple results dataframe
triple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(threeway)) {
  marker1<-threeway[[i]][1]
  marker2<-threeway[[i]][2]
  marker3<-threeway[[i]][3]
  test<-choose_three(genus_df, marker1, marker2, marker3)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  triple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, uniques=length(unique(colnames(y))))
  triple_results<-rbind(triple_results,triple_results.tmp)
}

# sorted by descending number of uniques
triple_sorted_genus<-triple_results[order(-triple_results$uniques),]

#reformat then merge with genus
triple_sorted_genus$marker<-paste(triple_sorted_genus$marker1, triple_sorted_genus$marker2, triple_sorted_genus$marker3)
triple_sorted_genus<-triple_sorted_genus[,-c(1:3)]
triple_sorted_genus$combo<-"3"
genus<-rbind(genus, triple_sorted_genus)

# initialize quadruple results dataframe
quadruple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fourway)) {
  marker1<-fourway[[i]][1]
  marker2<-fourway[[i]][2]
  marker3<-fourway[[i]][3]
  marker4<-fourway[[i]][4]
  test<-choose_four(genus_df, marker1, marker2, marker3, marker4)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quadruple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, uniques=length(unique(colnames(y))))
  quadruple_results<-rbind(quadruple_results,quadruple_results.tmp)
}

# sorted by descending number of uniques
quadruple_sorted_genus<-quadruple_results[order(-quadruple_results$uniques),]

#reformat then merge with genus
quadruple_sorted_genus$marker<-paste(quadruple_sorted_genus$marker1, quadruple_sorted_genus$marker2, quadruple_sorted_genus$marker3, quadruple_sorted_genus$marker4)
quadruple_sorted_genus<-quadruple_sorted_genus[,-c(1:4)]
quadruple_sorted_genus$combo<-"4"
genus<-rbind(genus, quadruple_sorted_genus)

# initialize quintuple results dataframe
quintuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fiveway)) {
  marker1<-fiveway[[i]][1]
  marker2<-fiveway[[i]][2]
  marker3<-fiveway[[i]][3]
  marker4<-fiveway[[i]][4]
  marker5<-fiveway[[i]][5]
  
  test<-choose_five(genus_df, marker1, marker2, marker3, marker4, marker5)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quintuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, uniques=length(unique(colnames(y))))
  quintuple_results<-rbind(quintuple_results,quintuple_results.tmp)
}

# sorted by descending number of uniques
quintuple_sorted_genus<-quintuple_results[order(-quintuple_results$uniques),]

#reformat then merge with genus
quintuple_sorted_genus$marker<-paste(quintuple_sorted_genus$marker1, quintuple_sorted_genus$marker2, quintuple_sorted_genus$marker3, quintuple_sorted_genus$marker4, quintuple_sorted_genus$marker5)
quintuple_sorted_genus<-quintuple_sorted_genus[,-c(1:5)]
quintuple_sorted_genus$combo<-"5"
genus<-rbind(genus, quintuple_sorted_genus)

# initialize sextuple results dataframe
sextuple_results<-data.frame(marker1=character(),
                             marker2=character(),
                             marker3=character(),
                             marker4=character(),
                             marker5=character(),
                             marker6=character(),
                             uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(sixway)) {
  marker1<-sixway[[i]][1]
  marker2<-sixway[[i]][2]
  marker3<-sixway[[i]][3]
  marker4<-sixway[[i]][4]
  marker5<-sixway[[i]][5]
  marker6<-sixway[[i]][6]
  
  test<-choose_six(genus_df, marker1, marker2, marker3, marker4, marker5, marker6)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  sextuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, marker6=marker6, uniques=length(unique(colnames(y))))
  sextuple_results<-rbind(sextuple_results, sextuple_results.tmp)
}

# sorted by descending number of uniques
sextuple_sorted_genus<-sextuple_results[order(-sextuple_results$uniques),]

#reformat then merge with genus
sextuple_sorted_genus$marker<-paste(sextuple_sorted_genus$marker1, sextuple_sorted_genus$marker2, sextuple_sorted_genus$marker3, sextuple_sorted_genus$marker4, sextuple_sorted_genus$marker5, sextuple_sorted_genus$marker6)
sextuple_sorted_genus<-sextuple_sorted_genus[,-c(1:6)]
sextuple_sorted_genus$combo<-"6"
genus<-rbind(genus, sextuple_sorted_genus)
#create copy for editing
genus2<-genus

# replace letter codes with primer names
genus2$marker<-gsub("F", "fwh1", genus2$marker)
genus2$marker<-gsub("B", "F230R", genus2$marker)
genus2$marker<-gsub("A","BR5",genus2$marker)
genus2$marker<-gsub("C", "ml-jg", genus2$marker)
genus2$marker<-gsub("D", "BF1R2", genus2$marker)
genus2$marker<-gsub("E", "BF2R2", genus2$marker)

# Plot an accumulation curve
f3<-ggplot(genus2, aes(x=combo, y=uniques, label = marker)) +
  ggtitle("3)") +
  geom_point() +
  annotate("text", label="A\n", x=1, y=55, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="B + C\n", x=2, y=70, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C", x=3, y=85, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + F", x=4, y=90, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ F", x=5, y=90, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E + F", x=6, y=90, hjust=0.5, vjust=-0.5, size=2.25) +
  labs(x="No. Combined Markers", y="Genus richness") +
  ylim(1,125) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

##############################################
# Create tables for family rank
##############################################

# initialize single primer results dataframe
single_results<-data.frame(marker=character(),
                           uniques=integer())

# get number of unique taxa per marker
for (i in markers) {
  p1<-paste("^",i,"_",sep="")
  target<-rownames(family_df)
  x<-family_df[grep(p1,target),]
  y<-x[,colSums(x) !=0]
  single_results.tmp<-data.frame(marker=i, uniques=length(unique(colnames(y))))
  single_results<-rbind(single_results,single_results.tmp)
}

# sorted by descending number of uniques
single_sorted_family<-single_results[order(-single_results$uniques),]

# add column to indicate number of markers combined
single_sorted_family$combo<-"1"

# initialize double results dataframe
double_results<-data.frame(marker1=character(),
                           marker2=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(pairwise)) {
  marker1<-pairwise[[i]][1]
  marker2<-pairwise[[i]][2]
  test<-choose_two(family_df, marker1, marker2)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  double_results.tmp<-data.frame(marker1=marker1, marker2=marker2, uniques=length(unique(colnames(y))))
  double_results<-rbind(double_results,double_results.tmp)
}

# sorted by descending number of uniques
double_sorted_family<-double_results[order(-double_results$uniques),]

#reformat then merge with single_sorted_family
double_sorted_family$marker<-paste(double_sorted_family$marker1, double_sorted_family$marker2)
double_sorted_family<-double_sorted_family[,-c(1:2)]
double_sorted_family$combo<-"2"
family<-rbind(single_sorted_family, double_sorted_family)

# initialize triple results dataframe
triple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(threeway)) {
  marker1<-threeway[[i]][1]
  marker2<-threeway[[i]][2]
  marker3<-threeway[[i]][3]
  test<-choose_three(family_df, marker1, marker2, marker3)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  triple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, uniques=length(unique(colnames(y))))
  triple_results<-rbind(triple_results,triple_results.tmp)
}

# sorted by descending number of uniques
triple_sorted_family<-triple_results[order(-triple_results$uniques),]

#reformat then merge with family
triple_sorted_family$marker<-paste(triple_sorted_family$marker1, triple_sorted_family$marker2, triple_sorted_family$marker3)
triple_sorted_family<-triple_sorted_family[,-c(1:3)]
triple_sorted_family$combo<-"3"
family<-rbind(family, triple_sorted_family)

# initialize quadruple results dataframe
quadruple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fourway)) {
  marker1<-fourway[[i]][1]
  marker2<-fourway[[i]][2]
  marker3<-fourway[[i]][3]
  marker4<-fourway[[i]][4]
  test<-choose_four(family_df, marker1, marker2, marker3, marker4)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quadruple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, uniques=length(unique(colnames(y))))
  quadruple_results<-rbind(quadruple_results,quadruple_results.tmp)
}

# sorted by descending number of uniques
quadruple_sorted_family<-quadruple_results[order(-quadruple_results$uniques),]

#reformat then merge with family
quadruple_sorted_family$marker<-paste(quadruple_sorted_family$marker1, quadruple_sorted_family$marker2, quadruple_sorted_family$marker3, quadruple_sorted_family$marker4)
quadruple_sorted_family<-quadruple_sorted_family[,-c(1:4)]
quadruple_sorted_family$combo<-"4"
family<-rbind(family, quadruple_sorted_family)

# initialize quintuple results dataframe
quintuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fiveway)) {
  marker1<-fiveway[[i]][1]
  marker2<-fiveway[[i]][2]
  marker3<-fiveway[[i]][3]
  marker4<-fiveway[[i]][4]
  marker5<-fiveway[[i]][5]
  
  test<-choose_five(family_df, marker1, marker2, marker3, marker4, marker5)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quintuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, uniques=length(unique(colnames(y))))
  quintuple_results<-rbind(quintuple_results,quintuple_results.tmp)
}

# sorted by descending number of uniques
quintuple_sorted_family<-quintuple_results[order(-quintuple_results$uniques),]

#reformat then merge with family
quintuple_sorted_family$marker<-paste(quintuple_sorted_family$marker1, quintuple_sorted_family$marker2, quintuple_sorted_family$marker3, quintuple_sorted_family$marker4, quintuple_sorted_family$marker5)
quintuple_sorted_family<-quintuple_sorted_family[,-c(1:5)]
quintuple_sorted_family$combo<-"5"
family<-rbind(family, quintuple_sorted_family)

# initialize sextuple results dataframe
sextuple_results<-data.frame(marker1=character(),
                             marker2=character(),
                             marker3=character(),
                             marker4=character(),
                             marker5=character(),
                             marker6=character(),
                             uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(sixway)) {
  marker1<-sixway[[i]][1]
  marker2<-sixway[[i]][2]
  marker3<-sixway[[i]][3]
  marker4<-sixway[[i]][4]
  marker5<-sixway[[i]][5]
  marker6<-sixway[[i]][6]
  
  test<-choose_six(family_df, marker1, marker2, marker3, marker4, marker5, marker6)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  sextuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, marker6=marker6, uniques=length(unique(colnames(y))))
  sextuple_results<-rbind(sextuple_results, sextuple_results.tmp)
}

# sorted by descending number of uniques
sextuple_sorted_family<-sextuple_results[order(-sextuple_results$uniques),]

#reformat then merge with family
sextuple_sorted_family$marker<-paste(sextuple_sorted_family$marker1, sextuple_sorted_family$marker2, sextuple_sorted_family$marker3, sextuple_sorted_family$marker4, sextuple_sorted_family$marker5, sextuple_sorted_family$marker6)
sextuple_sorted_family<-sextuple_sorted_family[,-c(1:6)]
sextuple_sorted_family$combo<-"6"
family<-rbind(family, sextuple_sorted_family)
#create copy for editing
family2<-family

# replace letter codes with primer names
family2$marker<-gsub("F", "fwh1", family2$marker)
family2$marker<-gsub("B", "F230R", family2$marker)
family2$marker<-gsub("A","BR5",family2$marker)
family2$marker<-gsub("C", "ml-jg", family2$marker)
family2$marker<-gsub("D", "BF1R2", family2$marker)
family2$marker<-gsub("E", "BF2R2", family2$marker)

# Plot an accumulation curve
f4<-ggplot(family2, aes(x=combo, y=uniques, label = marker)) +
  ggtitle("4)") +
  geom_point() +
  annotate("text", label="A / B\n", x=1, y=30, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B /\nB + C", x=2, y=45, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="B + C + E /\nB + C + F", x=3, y=50, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="B + C + E + F", x=4, y=57, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + E\n+ F", x=5, y=57, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D\n+ E + F", x=6, y=58, hjust=0.5, vjust=-0.5, size=2.25) +
  labs(x="No. Combined Markers", y="Family richness") +
  ylim(1,100) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

##############################################
# Create tables for order rank
##############################################

# initialize single primer results dataframe
single_results<-data.frame(marker=character(),
                           uniques=integer())

# get number of unique taxa per marker
for (i in markers) {
  p1<-paste("^",i,"_",sep="")
  target<-rownames(order_df)
  x<-order_df[grep(p1,target),]
  y<-x[,colSums(x) !=0]
  single_results.tmp<-data.frame(marker=i, uniques=length(unique(colnames(y))))
  single_results<-rbind(single_results,single_results.tmp)
}

# sorted by descending number of uniques
single_sorted_order<-single_results[order(-single_results$uniques),]

# add column to indicate number of markers combined
single_sorted_order$combo<-"1"

# initialize double results dataframe
double_results<-data.frame(marker1=character(),
                           marker2=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(pairwise)) {
  marker1<-pairwise[[i]][1]
  marker2<-pairwise[[i]][2]
  test<-choose_two(order_df, marker1, marker2)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  double_results.tmp<-data.frame(marker1=marker1, marker2=marker2, uniques=length(unique(colnames(y))))
  double_results<-rbind(double_results,double_results.tmp)
}

# sorted by descending number of uniques
double_sorted_order<-double_results[order(-double_results$uniques),]

#reformat then merge with single_sorted_order
double_sorted_order$marker<-paste(double_sorted_order$marker1, double_sorted_order$marker2)
double_sorted_order<-double_sorted_order[,-c(1:2)]
double_sorted_order$combo<-"2"
order<-rbind(single_sorted_order, double_sorted_order)

# initialize triple results dataframe
triple_results<-data.frame(marker1=character(),
                           marker2=character(),
                           marker3=character(),
                           uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(threeway)) {
  marker1<-threeway[[i]][1]
  marker2<-threeway[[i]][2]
  marker3<-threeway[[i]][3]
  test<-choose_three(order_df, marker1, marker2, marker3)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  triple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, uniques=length(unique(colnames(y))))
  triple_results<-rbind(triple_results,triple_results.tmp)
}

# sorted by descending number of uniques
triple_sorted_order<-triple_results[order(-triple_results$uniques),]

#reformat then merge with order
triple_sorted_order$marker<-paste(triple_sorted_order$marker1, triple_sorted_order$marker2, triple_sorted_order$marker3)
triple_sorted_order<-triple_sorted_order[,-c(1:3)]
triple_sorted_order$combo<-"3"
order<-rbind(order, triple_sorted_order)

# initialize quadruple results dataframe
quadruple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fourway)) {
  marker1<-fourway[[i]][1]
  marker2<-fourway[[i]][2]
  marker3<-fourway[[i]][3]
  marker4<-fourway[[i]][4]
  test<-choose_four(order_df, marker1, marker2, marker3, marker4)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quadruple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, uniques=length(unique(colnames(y))))
  quadruple_results<-rbind(quadruple_results,quadruple_results.tmp)
}

# sorted by descending number of uniques
quadruple_sorted_order<-quadruple_results[order(-quadruple_results$uniques),]

#reformat then merge with order
quadruple_sorted_order$marker<-paste(quadruple_sorted_order$marker1, quadruple_sorted_order$marker2, quadruple_sorted_order$marker3, quadruple_sorted_order$marker4)
quadruple_sorted_order<-quadruple_sorted_order[,-c(1:4)]
quadruple_sorted_order$combo<-"4"
order<-rbind(order, quadruple_sorted_order)

# initialize quintuple results dataframe
quintuple_results<-data.frame(marker1=character(),
                              marker2=character(),
                              marker3=character(),
                              marker4=character(),
                              marker5=character(),
                              uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(fiveway)) {
  marker1<-fiveway[[i]][1]
  marker2<-fiveway[[i]][2]
  marker3<-fiveway[[i]][3]
  marker4<-fiveway[[i]][4]
  marker5<-fiveway[[i]][5]
  
  test<-choose_five(order_df, marker1, marker2, marker3, marker4, marker5)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  quintuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, uniques=length(unique(colnames(y))))
  quintuple_results<-rbind(quintuple_results,quintuple_results.tmp)
}

# sorted by descending number of uniques
quintuple_sorted_order<-quintuple_results[order(-quintuple_results$uniques),]

#reformat then merge with order
quintuple_sorted_order$marker<-paste(quintuple_sorted_order$marker1, quintuple_sorted_order$marker2, quintuple_sorted_order$marker3, quintuple_sorted_order$marker4, quintuple_sorted_order$marker5)
quintuple_sorted_order<-quintuple_sorted_order[,-c(1:5)]
quintuple_sorted_order$combo<-"5"
order<-rbind(order, quintuple_sorted_order)

# initialize sextuple results dataframe
sextuple_results<-data.frame(marker1=character(),
                             marker2=character(),
                             marker3=character(),
                             marker4=character(),
                             marker5=character(),
                             marker6=character(),
                             uniques=integer())

## for each column in pairwise, run combo_uniques
for (i in names(sixway)) {
  marker1<-sixway[[i]][1]
  marker2<-sixway[[i]][2]
  marker3<-sixway[[i]][3]
  marker4<-sixway[[i]][4]
  marker5<-sixway[[i]][5]
  marker6<-sixway[[i]][6]
  
  test<-choose_six(order_df, marker1, marker2, marker3, marker4, marker5, marker6)
  x<-test$matrix
  y<-x[,colSums(x) !=0]
  sextuple_results.tmp<-data.frame(marker1=marker1, marker2=marker2, marker3=marker3, marker4=marker4, marker5=marker5, marker6=marker6, uniques=length(unique(colnames(y))))
  sextuple_results<-rbind(sextuple_results, sextuple_results.tmp)
}

# sorted by descending number of uniques
sextuple_sorted_order<-sextuple_results[order(-sextuple_results$uniques),]

#reformat then merge with order
sextuple_sorted_order$marker<-paste(sextuple_sorted_order$marker1, sextuple_sorted_order$marker2, sextuple_sorted_order$marker3, sextuple_sorted_order$marker4, sextuple_sorted_order$marker5, sextuple_sorted_order$marker6)
sextuple_sorted_order<-sextuple_sorted_order[,-c(1:6)]
sextuple_sorted_order$combo<-"6"
order<-rbind(order, sextuple_sorted_order)
#create copy for editing
order2<-order

# replace letter codes with primer names
order2$marker<-gsub("F", "fwh1", order2$marker)
order2$marker<-gsub("B", "F230R", order2$marker)
order2$marker<-gsub("A","BR5",order2$marker)
order2$marker<-gsub("C", "ml-jg", order2$marker)
order2$marker<-gsub("D", "BF1R2", order2$marker)
order2$marker<-gsub("E", "BF2R2", order2$marker)

# Plot an accumulation curve
f5<-ggplot(order2, aes(x=combo, y=uniques, label = marker)) +
  ggtitle("5)") +
  geom_point() +
  annotate("text", label="D\n", x=1, y=27, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B /\nA + D /\nB + D /\nD + F", x=2, y=30, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + D + F", x=3, y=37, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + D + F /\nA + C + D + F", x=4, y=39, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D \n+ F", x=5, y=40, hjust=0.5, vjust=-0.5, size=2.25) +
  annotate("text", label="A + B + C + D \n+ E + F", x=6, y=40, hjust=0.5, vjust=-0.5, size=2.25) +
  labs(x="No. Combined Markers", y="Order richness") +
  ylim(1,50) + 
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))

# print both plots side by side
p<-grid.arrange(f1, f2, f3, f4, f5, ncol=2)

ggsave("F2_multimarker_richness_accumulation.pdf", p, width=8.5, height=11, units="in")
