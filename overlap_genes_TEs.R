library(data.table)
library(stats4)
library(BiocGenerics)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)
library(GenomicRanges)
library(intervals)
library(purrr)

#Import tables

Genes=read.csv("abc.csv", header = F, sep = "\t")
head(Genes)

TEs=read.csv("dari_03-02_filtered_size_merged_insertions.out", header = T ,sep = ",")
names(TEs) <- c("V1", "V1", "V3", "V4", "V5", "Scaffold","Start","End","Length","V10","TE_subfamily","TE_class")
head(TEs)
dim(TEs)

TEs_spec = TEs[, c("Scaffold", "Start", "End", "Subfamily")]
head(TEs_spec)

##### STEP 1 #####
#Create intervals to genes and TEs
GeneRanges=GRanges(seqnames=Genes$V1,ranges=IRanges(start=as.numeric(as.vector(Genes$V4)),end=as.numeric(as.vector(Genes$V5))), name=Genes$V9)
#TEsRanges=GRanges(seqnames=TEs$Scaffold,ranges=IRanges(start=as.numeric(as.vector(TEs$Start)),end=as.numeric(as.vector(TEs$End))),name=TEs$Class.family)

TEsRanges=GRanges(seqnames=TEs_spec$Scaffold,ranges=IRanges(start=as.numeric(as.vector(TEs_spec$Start)),end=as.numeric(as.vector(TEs_spec$End))),name=TEs_spec$Subfamily)

#Search for overlap sequences 
Overlaps=findOverlapPairs(GeneRanges, TEsRanges)
#PosTable=Overlaps@second # coordinates of the TEs that are inside the Genes
#PosTable=Overlaps@second$name # IDs of the ones that are inside

#Save overlapped sequences
write.table(Overlaps, "Overlaps.csv", sep = "\t")
#write.table(PosTable, "TEs_inside_genes.csv")


#Simplest table

Overlaps = read.csv("Overlaps.csv", header = T, sep = "\t")
head(Overlaps)
dim(Overlaps)
IDs_genes_TEs = data.frame(Overlaps$first.name, Overlaps$first.X.width ,Overlaps$second.nam, Overlaps$second.X.width, Overlaps$second.X.seqnames ,Overlaps$second.X.start, Overlaps$second.X.end)
names(IDs_genes_TEs) <- c("Gene_ID", "Gene_length", "TE_subfamily", "TE_length", "Scaffold", "Start", "End")
write.table(IDs_genes_TEs, "Genes_with_TEs_inside.csv", sep ="\t",row.names = F )

##### STEP 2 #####

# To do the Up and Downstream what I did was just change the start and end of the gene + or - the quantity you want
x = 3000
Genes$Start_UP=Genes$V4-x
Genes$End_Down=Genes$V5+x

head(Genes) # Are the columns there?


#Create range for upstream gene position
GenesRanges_UP=GRanges(seqnames=Genes$V1,ranges=IRanges(start=as.numeric(as.vector(Genes$Start_UP)),end=as.numeric(as.vector(Genes$V4))),name=Genes$V9)

#Create range for downstream gene position
GenesRanges_DOWN=GRanges(seqnames=Genes$V1,ranges=IRanges(start=as.numeric(as.vector(Genes$V5)),end=as.numeric(as.vector(Genes$End_Down))),name=Genes$V9)

##### STEP 3 #####

#Save table with TEs sequences at upstream region
TEs_upstream = findOverlapPairs(GenesRanges_UP, TEsRanges)
write.table(TEs_upstream, "Genes_with_TEs_2kb_UPstream.csv", sep = "\t")                          
                          
#Save table with TEs sequences at downstream region
TEs_downstream = findOverlapPairs(GenesRanges_DOWN, TEsRanges)
write.table(TEs_downstream, "Genes_with_TEs_2kb_DOWNstream.csv", sep = "\t")

#Simplest table


TE_up = read.csv("Genes_with_TEs_2kb_UPstream.csv", header = T, sep = "\t")
TE_up_df = data.frame(TE_up$first.name, TE_up$first.X.width ,TE_up$second.name, TE_up$second.X.width, TE_up$second.X.seqnames ,TE_up$second.X.start, TE_up$second.X.end)
head(TE_up_df)
names(TE_up_df) <- c("Gene ID", "TE_subfamily","TE length","Scaffold", "Start", "End")
write.table(TE_up_df, "Genes_with_TEs_2kb_UPstream.csv", sep ="\t",row.names = F)

TE_down = read.csv("Genes_with_TEs_2kb_DOWNstream.csv", header = T, sep = "\t")
TE_down_df = data.frame(TE_down$first.name, TE_down$first.X.width ,TE_down$second.name, TE_down$second.X.width, TE_down$second.X.seqnames ,TE_down$second.X.start, TE_down$second.X.end)
names(TE_down_df) <- c("Gene ID", "TE_subfamily", "TE length","Scaffold", "Start", "End")
head(TE_down_df)
write.table(TE_down_df, "Genes_with_TEs_2kb_DOWNstream.csv", sep ="\t",row.names = F)

