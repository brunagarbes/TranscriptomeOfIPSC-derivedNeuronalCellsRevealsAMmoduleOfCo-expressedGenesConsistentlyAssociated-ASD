###NormalizationStepsNPCs###

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);

#load required packages
library(RUVSeq)
library(EDASeq)
library(WGCNA)
library(calibrate)
library(DESeq2)
library(pheatmap)
library (RColorBrewer)
#You might have to download additional packages depending on your currently R library

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#loading input data

geneRPKM=read.delim ("allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneRPKM)=geneRPKM[,1]
geneRPKM=geneRPKM[,-1]

#load countdata
counts=read.table("countdata_20M_NPC.txt", header=TRUE)
sampleInfo=read.delim("samplesheet_NPC.txt", sep="\t", header=TRUE)
HKgenes_full <- read.csv ("HK_full_gene_list_ensembl_biomart.csv")

#removing samples with >55%
#samples with neuronal cell proportion > 55% are removed
countData <- counts [,c(1:10,12:14,22:37)]
sampleData <- sampleInfo [c(1:10,12:14,22:37),]
sampleData <- sampleData[,c(1:2,4,9:10)]
geneRPKM2 <- geneRPKM [,c(1:10,12:14,22:37)]

filtered <- which(rowSums(geneRPKM2 > 1) >= 15)
RPKM <- geneRPKM2[filtered,]
countData=subset (countData, rownames(countData) %in% rownames(RPKM))


#find the intersect between countData and HKgenes
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData))


#first differential expression analysis with DESeq2
#we generated a gene list which is known to be no differentially expressed (pvalue>=0.7)
#so we can built up a housekeeping gene list for RUVseq analysis
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)

#set which group is the reference group
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
#you can check the number of outliers and lowcount genes using the function summary
summary(res)
res <- as.data.frame(res)
nDEGs_deseq <- subset (res, res$padj>=0.7)
nDEGs_deseq <- rownames(nDEGs_deseq)

nDEGs_HKgenes <- intersect(nDEGs_deseq,HKgenes)
controlGenes <- nDEGs_HKgenes
m=match(colnames(countData), sampleData$label)

#RUVSeq

#before normalization
#PCA and hierarchical cluster graphs before normalization

pdf("BeforeNorm_plots_NPCs55.pdf")
boxplot(log2(countData))
pca=princomp(cor(countData, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,beforeNorm55, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("Neur_prop_group"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

pdf("contrGenes_nDEGs_HKgenes_NPCs55.pdf")

 for( i in 1:4){
  
  RUV <- RUVg(as.matrix(countData), controlGenes, k=i)
  N <- RUV$normalizedCounts
  
  boxplot(log2(N))
  title("Boxplot NPCs Samples, 1<= K <=4")
  pca=princomp(cor(N, method="s"))
  plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="RUV Normalization PCA graphs considering neuron proportions", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
  dds <- DESeqDataSetFromMatrix(countData = N, colData = sampleData, design = ~ condition)
  vsd <- varianceStabilizingTransformation(dds)
  plotPCA(vsd, intgroup=c("Neur_prop_group"))
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(colnames(countData))
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
 }
dev.off()

###DESeq2 - Differential expression analysis###
#create deseq object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition) #so preciso criar o dds com o countdata

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)

###WGCNA analysis###
#check for variance over the mean and possibliy remove some genes
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)

keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.020)
#we calculated the variance over the mean and tried different cuttoffs to see how many genes were still retained
ncvsd=ncvsd[keep,]

#transpose the table
dat=t(ncvsd)
infoData=rownames(ncvsd)

#run WGCNA#
#run soft pick threshold
#test a series of powers to which co-expression similarity is raised to calculate adjacency
powers=c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, blockSize = 20000, networkType = "signed") # Signed gives an indication of positive vs negative correlations
#plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

pdf(file="WGCNA_SoftThreshold_NPCs_nDEGs_HKgenes_p0_7_k2.pdf")
#fit to scale-free topology
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red");
abline(h=c(0.8, 0.5), col="red")

#mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col="red")
dev.off()
