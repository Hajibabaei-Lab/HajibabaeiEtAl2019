# Teresita M. Porter, May 30, 2019

library(vegan)
library(reshape2)
library(gridExtra)
library(grid)
library(ggplot2)
library(plyr)
library(scales) # view default colours

###################################################################

# Read infile prepared by python script
A<-read.table(file="taxonomy.matrix", head=TRUE, sep=",")

# Summarize ESVs in all detected phyla
A2 <- dcast(A, Phylum ~ . , value.var="A_GlobalESV", fun.aggregate=length)
names(A2)<-c("Phyla","ESVs")

# Sort descending
A2<-A2[order(-A2$ESVs),]

# Keep top 10
A2<-A2[1:10,]

# Summarize reads in all detected phyla
A3 <- dcast(A, Phylum ~ . , value.var="ESVsize", fun.aggregate=sum)
names(A3)<-c("Phyla","Reads")

# Sort descending
A3<-A3[order(-A3$Reads),]

# Keep top 10
A3<-A3[1:10,]

# Create phylum summary table
phylumTable<-merge(A2,A3,by="Phyla")

# Create long form form ggplot
phylumTable.long<-melt(phylumTable, id=c("Phyla"))

# Create 100% stacked horizontal bar plot to summarize phyla diversity
p1<-ggplot(data=phylumTable.long, aes(x=variable, y=value, fill=Phyla)) +
  ggtitle("A) All phyla") +
  geom_bar(position="fill", stat="identity") +
  labs(y="Proportion") +
  scale_y_continuous(labels = percent_format()) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),
        axis.title.x=element_blank()) +
  # annotate top 3 ESVs
  annotate("text",x=1, y=0.8,label="11,285") +
  annotate("text",x=1, y=0.3,label="9,862") +
  annotate("text",x=1, y=0.07,label="5,012") +
  # annotate top 3 Reads
  annotate("text",x=2, y=0.55,label="1,280,397") +
  annotate("text",x=2, y=0.92,label="358,274") +
  annotate("text",x=2, y=0.08,label="237,614") 

# Select phylum Arthropoda only
Arth_df<-A[A$Phylum=="Arthropoda",]

# Plot orders EPTC and Other
A4 <- dcast(Arth_df, Order+Family+fBP ~ . , value.var="ESVsize", fun.aggregate=sum)
names(A4)<-c("Order","Family", "fBP", "Counts")

A4_E<-A4[A4$Order=="Ephemeroptera",c(1,4)]
names(A4_E)<-c("Taxon","Counts")
A4_P<-A4[A4$Order=="Plecoptera_Insecta",c(1,4)]
names(A4_P)<-c("Taxon","Counts")
A4_T<-A4[A4$Order=="Trichoptera",c(1,4)]
names(A4_T)<-c("Taxon","Counts")
A4_C<-A4[A4$Family=="Chironomidae" &
         A4$fBP>=0.20, c(2,4)]
names(A4_C)<-c("Taxon","Counts")

# remove target assemblage
A4b<-A4[!(A4$Order=="Ephemeroptera" |
         A4$Order=="Plecoptera_Insecta" |
         A4$Order=="Trichoptera" |
         A4$Family=="Chironomidae"),]

# Calc what's left 'other'
other<-sum(A4b$Counts)
A4_O<-data.frame(Taxon="Other",Counts=other)

# Put it back together
A5<-rbind(A4_E,A4_P,A4_T,A4_C,A4_O)
A5$rank<-"Reads"

A6 <- dcast(Arth_df, Order+Family+fBP ~ . , value.var="A_GlobalESV", fun.aggregate = function(x){length(unique(x))})
names(A6)<-c("Order","Family","fBP", "Counts")

A6_E<-A6[A6$Order=="Ephemeroptera",c(1,4)]
names(A6_E)<-c("Taxon","Counts")
A6_P<-A6[A6$Order=="Plecoptera_Insecta",c(1,4)]
names(A6_P)<-c("Taxon","Counts")
A6_T<-A6[A6$Order=="Trichoptera",c(1,4)]
names(A6_T)<-c("Taxon","Counts")
A6_C<-A6[A6$Family=="Chironomidae" &
         A6$fBP>=0.20, c(2,4)]
names(A6_C)<-c("Taxon","Counts")

# Remove target assemblage
A6b<-A6[!(A6$Order=="Ephemeroptera" |
          A6$Order=="Plecoptera_Insecta" |
          A6$Order=="Trichoptera" |
          A6$Family=="Chironomidae"),]

other<-sum(A6b$Counts)
A6_O<-data.frame(Taxon="Other",Counts=other)

A7<-rbind(A6_E,A6_P,A6_T,A6_C,A6_O)
A7$rank<-"ESVs"

orderTable.long<-rbind(A5,A7)

# view colour codes
show_col(hue_pal()(4))

# Create 100% stacked horizontal bar plot to summarize order diversity focusing on EPTC
p2<-ggplot(data=orderTable.long, aes(x=rank, y=Counts, fill=Taxon)) +
  ggtitle("B) Arthropoda taxa") +
  geom_bar(position="fill", stat="identity") +
  labs(y="Proportion", x="Taxa") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values=c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF", "grey")) +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        text=element_text(size=14),
        axis.title.x=element_blank()) 

g<-grid.arrange(p1,p2,ncol=1)
ggsave("FigS2_TaxSummary.pdf",g, width=8.5, height=8.5, units="in")