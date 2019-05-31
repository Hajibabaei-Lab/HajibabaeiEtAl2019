# Teresita M. Porter, May 30, 2019

library(stringr)
library(scales)
library(reshape2)
library(vegan)
library(RColorBrewer)

###################################################################
# Edit rarecurve function to remove the horizontal lines
###################################################################

rarecurve2 <- function (x, step = 1, sample, xlab = "Sample Size", ylab = "Species", 
                        label = TRUE, col, lty, ...) 
{
  x <- as.matrix(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if (missing(col)) 
    col <- par("col")
  if (missing(lty)) 
    lty <- par("lty")
  tot <- rowSums(x)
  S <- specnumber(x)
  if (any(S <= 0)) {
    message("empty rows removed")
    x <- x[S > 0, , drop = FALSE]
    tot <- tot[S > 0]
    S <- S[S > 0]
  }
  nr <- nrow(x)
  col <- rep(col, length.out = nr)
  lty <- rep(lty, length.out = nr)
  out <- lapply(seq_len(nr), function(i) {
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) 
      n <- c(n, tot[i])
    drop(rarefy(x[i, ], n))
  })
  Nmax <- sapply(out, function(x) max(attr(x, "Subsample")))
  Smax <- sapply(out, max)
  plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = xlab, ylab = ylab, 
       type = "n", ...)
  if (!missing(sample)) {
    abline(v = sample) 
    rare <- sapply(out, function(z) approx(x = attr(z, "Subsample"), 
                                           y = z, xout = sample, rule = 1)$y)
    #    abline(h = rare, lwd = 0.5) #turn off horizontal lines
  }
  for (ln in seq_along(out)) {
    N <- attr(out[[ln]], "Subsample")
    lines(N, out[[ln]], col = col[ln], lty = lty[ln], ...)
  }
  if (label) { 
    ordilabel(cbind(tot, S), labels = rownames(x), ...)
  }
  invisible(out)
}

#####################################################################

# Read in sample x taxonomy table
A<-read.csv(file="taxonomy.matrix", head=TRUE)

# Filter table for Arthropoda only
B<-A[A$Phylum=="Arthropoda",]

# Filter table for species with bp >= 0.70 for 95% correct
C<-B[B$sBP>=0.70,]

# Filter table for genera with bp >= 0.3
D<-B[B$gBP>=0.30,]

# Filter table for families with bp >= 0.2
E<-B[B$fBP>=0.30,]

###
# Figureout proportion of Arthropoda species that are confidently identified
###

# all arthropoda species
x<-length(unique(B$Species[grep("SoilKit", B$SampleName)]))
# all arthropoda species with sBP>=0.70
y<-length(unique(C$Species[grep("SoilKit", C$SampleName)]))

# all arthropoda genera
x2<-length(unique(B$Genus[grep("SoilKit", B$SampleName)]))
# all arthropoda genera with sBP>=0,30
y2<-length(unique(D$Genus[grep("SoilKit", D$SampleName)]))

# all arthropoda families
x3<-length(unique(B$Family[grep("SoilKit", B$SampleName)]))
# all arthropoda genera with sBP>=0,30
y3<-length(unique(E$Genus[grep("SoilKit", E$SampleName)]))


#######################################################
# Create dataframe for rarefaction curves
## pool data across reps
## separate curves for primers and kits and sites
## plot at esv, species, genus, and family ranks
######################################################

# Split up SampleName for Arthropoda matrix (for esv, order, class ranks)
B2<-data.frame(B, do.call(rbind, str_split(B$SampleName,"_")))

# Filter for SoilKit only
B2<-B2[grepl("SoilKit",B2$SampleName),]

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

# pivot to make esv, order, class matrices
esv<-dcast(B3, kit+primer+site+rep ~ A_GlobalESV, value.var = "ESVsize", fun.aggregate = sum)
order<-dcast(B3, kit+primer+site+rep ~ Order, value.var = "ESVsize", fun.aggregate = sum)
class<-dcast(B3, kit+primer+site+rep ~ Class, value.var = "ESVsize", fun.aggregate = sum)
# pivot to make species matrix sBP>=0.50 for 95% correct
B3_species<-B3[B3$sBP>=0.70,]
species<-dcast(B3_species, kit+primer+site+rep ~ Species, value.var = "ESVsize", fun.aggregate = sum)
# pivot to make genus matrix gBP>=0.30
B3_genus<-B3[B3$gBP>=0.30,]
genus<-dcast(B3_genus, kit+primer+site+rep ~ Genus, value.var = "ESVsize", fun.aggregate = sum)
# pivot to make family matrix fBP>=0.20
B3_family<-B3[B3$fBP>=0.20,]
family<-dcast(B3_family, kit+primer+site+rep ~ Family, value.var = "ESVsize", fun.aggregate = sum)

# merge kit+site+rep+primer into single column
esv$sample<-paste(esv$kit, esv$primer, esv$site, esv$rep, sep="_")
species$sample<-paste(species$kit, species$primer, species$site, species$rep, sep="_")
genus$sample<-paste(genus$kit, genus$primer, genus$site, genus$rep, sep="_")
family$sample<-paste(family$kit, family$primer, family$site, family$rep, sep="_")
order$sample<-paste(order$kit, order$primer, order$site, order$rep, sep="_")
class$sample<-paste(class$kit, class$primer, class$site, class$rep, sep="_")

# delete duplicate columns
esv$kit<-NULL
esv$primer<-NULL
esv$site<-NULL
esv$rep<-NULL
species$kit<-NULL
species$primer<-NULL
species$site<-NULL
species$rep<-NULL
genus$kit<-NULL
genus$primer<-NULL
genus$site<-NULL
genus$rep<-NULL
family$kit<-NULL
family$primer<-NULL
family$site<-NULL
family$rep<-NULL
order$kit<-NULL
order$primer<-NULL
order$site<-NULL
order$rep<-NULL
class$kit<-NULL
class$primer<-NULL
class$site<-NULL
class$rep<-NULL

#move SampleName to rownames then delete
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

#remove columns with only zeros
esv_notnull<-esv[,colSums(esv) !=0]
species_notnull<-species[,colSums(species) !=0]
genus_notnull<-genus[,colSums(genus) !=0]
family_notnull<-family[,colSums(family) != 0]
order_notnull<-order[,colSums(order) !=0]
class_notnull<-class[,colSums(class) !=0]

#remove rows with only zeros & edit rownames
esv_notnull2<-esv_notnull[rowSums(esv_notnull) !=0,]
species_notnull2<-species_notnull[rowSums(species_notnull) !=0,]
genus_notnull2<-genus_notnull[rowSums(genus_notnull) !=0,]
family_notnull2<-family_notnull[rowSums(family_notnull) !=0,]
order_notnull2<-order_notnull[rowSums(order_notnull) !=0,]
class_notnull2<-class_notnull[rowSums(class_notnull) !=0,]

#calculate 15th percentile for rrarefy function
esv_15percentile<-quantile(rowSums(esv_notnull2), prob=0.15)
species_15percentile<-quantile(rowSums(species_notnull2), prob=0.15)
genus_15percentile<-quantile(rowSums(genus_notnull2), prob=0.15)
family_15percentile<-quantile(rowSums(family_notnull2), prob=0.15)
order_15percentile<-quantile(rowSums(order_notnull2), prob=0.15)
class_15percentile<-quantile(rowSums(class_notnull2), prob=0.15)

###################################################################
##### Plot rarefaction curves
###################################################################

set.seed(1234)

#print rarefaction curves for each sample to assess read coverage per sample, indicate 15th percentile
pdf("FigS3_rarefaction.pdf")
par(mfrow=c(3,3),
    mar=c(4, 4, 3, 2))

esv_rarecurveout<-rarecurve2(esv_notnull2, sample=esv_15percentile, step=100, xlab="Reads", ylab="ESVs", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

species_rarecurveout<-rarecurve2(species_notnull2, sample=species_15percentile, step=100, xlab="Reads", ylab="Species", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

genus_rarecurveout<-rarecurve2(genus_notnull2, sample=genus_15percentile, step=100, xlab="Reads", ylab="Genus", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

family_rarecurveout<-rarecurve2(family_notnull2, sample=family_15percentile, step=100, xlab="Reads", ylab="Family", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

order_rarecurveout<-rarecurve2(order_notnull2, sample=order_15percentile, step=100, xlab="Reads", ylab="Order", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

class_rarecurveout<-rarecurve2(class_notnull2, sample=class_15percentile, step=100, xlab="Reads", ylab="Class", col=c(rep("#E41A1C",12),rep("#377EB8",12),rep("#4DAF4A",12),rep("#984EA3",12),rep("#FF7F00",12),rep("#A65628",12)), cex=0.6, cex.main=0.8, cex.lab=0.6, cex.axis=0.6, label=FALSE, cex.main=1.25, cex.lab=1.25, cex.axis=1.25, lty=c(rep(1,312)), lwd=1.5)

blank1<-plot(0,type='n',axes=FALSE,ann=FALSE)

blank2<-plot(0,type='n',axes=FALSE,ann=FALSE)

blank3<-plot(0,type='n',axes=FALSE,ann=FALSE)

blank3<-legend("topright", ncol=2, col=c("#E41A1C", "#377EB8","#4DAF4A","#984EA3","#FF7F00","#A65628"), lty=c(1,1,1,1,1,1), legend=c("BR5","F230R","ml-jg","BF1R2","BF2R2","fwh1"), bty="n", pt.cex=0.5, y.intersp=.75, lwd=1.5)

dev.off()

