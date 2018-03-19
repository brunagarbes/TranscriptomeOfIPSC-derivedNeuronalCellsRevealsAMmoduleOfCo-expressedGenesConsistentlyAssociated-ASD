
###NormalizationStepsNeurons###

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);

# Load required packages
library(RUVSeq)
library(EDASeq)
library(WGCNA)
library(calibrate)
library(DESeq2)
library(pheatmap)
library (RColorBrewer)
#You might have to download additional packages depending on your currently R library

#loading input data

geneFPKM=read.delim ("allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneFPKM)=geneFPKM[,1]
geneFPKM=geneFPKM[,c(41:63)]

counts=read.table("countdata_20M_neurons.txt", header=TRUE)
sampleData=read.delim("samplesheet_neurons.txt", sep="\t", header=TRUE)
HKgenes_full <- read.csv ("HK_full_gene_list_ensembl_biomart.csv")

#first filtering

filtered <- which(rowSums(geneFPKM > 1) >= 12)
FPKM <- geneFPKM[filtered,]
countData=subset (counts, rownames(counts) %in% rownames(FPKM))
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData))

#removing samples <50%
#samples with neuronal cell proportion <50% are removed

countData50 <- counts [,c(1,2,5:7,9:15,17:23)]
sampleData50 <- sampleData [c(1,2,5:7,9:15,17:23), c(1:5)]
geneFPKM50 <- geneFPKM [,c(1,2,5:7,9:15,17:23)]

filtered <- which(rowSums(geneFPKM50 > 1) >= 9)
FPKM <- geneFPKM50[filtered,]
countData50 =subset (countData50, rownames(countData50) %in% rownames(FPKM))

#first differential expression analysis
#differential expression with DESeq2
#we generated a gene list which is known to be no differentially expressed (pvalue>=0.4)
#so we can built up a housekeeping gene list for RUVseq analysis

dds <- DESeqDataSetFromMatrix(countData = countData50, colData = sampleData50, design = ~ condition)

#set which group is the reference group
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
#you can check the number of outliers and lowcount genes using the function summary
summary(res)
res <- as.data.frame(res)
nDEGs_deseq50 <- subset (res, res$pvalue>=0.4)
nDEGs_deseq50 <- rownames(nDEGs_deseq50)

#list of housekeeping genes used for the RUVseq normalization analysis
HKgenes50 <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData50))
nDEGs_HKgenes50 <- intersect(nDEGs_deseq50,HKgenes50)
m=match(colnames(countData50), sampleData50$label)


#RUV

#before normalization
#PCA and hierarchical cluster graphs before normalization

pdf("BeforeNorm_plots_neurons50.pdf")
boxplot(log2(countData50))
pca=princomp(cor(countData50, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,beforeNorm50, neur proportion", pch=20, col=labels2colors(sampleData50$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = countData50, colData = sampleData50, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("Neur_prop_group"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData50))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


#RUVseq analysis

controlGenes <- nDEGs_HKgenes50
pdf("contrGenes_nDEGs_HKgenes_neurons50.pdf")

for( i in 1:4){

 RUV <- RUVg(as.matrix(countData50), controlGenes, k=i)
 N <- RUV$normalizedCounts

 boxplot(log2(N))
 title("Boxplot Neuron Samples, 1<= K <=4")
 pca=princomp(cor(N, method="s"))
 plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="RUV Normalization PCA grph, neuron proportions", pch=20, col=labels2colors(sampleData50$Neur_prop_group[m]))
 dds <- DESeqDataSetFromMatrix(countData = N, colData = sampleData50, design = ~ condition)
 vsd <- varianceStabilizingTransformation(dds)
 plotPCA(vsd, intgroup=c("Neur_prop_group"))
 sampleDists <- dist(t(assay(vsd)))
 sampleDistMatrix <- as.matrix(sampleDists)
 rownames(sampleDistMatrix) <- paste(colnames(countData50))
 colnames(sampleDistMatrix) <- NULL
 colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
 pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
}
dev.off()

#RODAR WGCNA AQUI POIS ELE TAMBEM FOI USADO NA ESCOLHA DO MELHOR CONJUNTO DE DADOS
###DESeq2 - Differential expression analysis###
#create deseq object
dds <- DESeqDataSetFromMatrix(countData = countData50, colData = sampleData50, design = ~ condition) #so preciso criar o dds com o countdata

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)

###WGCNA analysis###
#check for variance over the mean and possibliy remove some genes
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)

keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.0125)
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
