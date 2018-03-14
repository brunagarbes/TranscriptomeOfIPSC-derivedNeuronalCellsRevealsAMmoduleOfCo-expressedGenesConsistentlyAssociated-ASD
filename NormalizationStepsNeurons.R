
###NormalizationStepsNeurons###



# Load required packages
library(RUVSeq)
library(EDASeq)
library(WGCNA)
library(calibrate)
library(DESeq2)
library(pheatmap)
library (RColorBrewer)

#load inputs/data
#load FPKM table to select genes in count table by FPKM>1
#load HKgenes lists
#load countdata and sample data

geneRPKM=read.delim ("D:/backup_arquivos_estagio/ruv-seq_analysis/allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneRPKM)=geneRPKM[,1]
geneRPKM=geneRPKM[,c(41:63)]

counts=read.table("D:/backup_arquivos_estagio/ruv-seq_analysis/countdata_20M_neurons.txt", header=TRUE)
sampleData=read.delim("D:/backup_arquivos_estagio/ruv-seq_analysis/samplesheet_neurons.txt", sep="\t", header=TRUE)
HKgenes_full <- read.csv ("D:/backup_arquivos_estagio/ruv-seq_analysis/HK_full_gene_list_ensembl_biomart.csv")


#filtering

filtered <- which(rowSums(geneRPKM > 1) >= 12)
RPKM <- geneRPKM[filtered,]
countData=subset (counts, rownames(counts) %in% rownames(RPKM))
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData))

#first past differential expression analysis

#differential expression with Deseq

dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
#set which group is the reference group
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
#you can check the number of outliers and lowcount genes using the function summary
summary(res)
res <- as.data.frame(res)
nDEGs_deseq <- subset (res, res$pvalue>=0.8)
nDEGs_deseq <- rownames(nDEGs_deseq)

#RUVseq

#before normalization

pdf("D:/backup_arquivos_estagio/ruv-seq_analysis/Results_RUVSeq_DROPBOX/BeforeNorm_plots_neurons.pdf")
boxplot(log2(countData))
pca=princomp(cor(countData, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,beforeNorm, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, intgroup=c("Neur_prop_group"))
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(countData))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap::pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()



#nDEG_HKgenes

#aqui o RUVseq � rodado com k=1 ate k=5 para testar a lista nDEGs_HKgenes que � uma combina��o dos genes n�o expressos das nossas amostras mais uma lista de genes knockout do paper de Eisenberg etal (2013)

controlGenes <- nDEGs_HKgenes
pdf("D:/backup_arquivos_estagio/ruv-seq_analysis/Results_RUVSeq_DROPBOX/contrGenes_nDEGs_HKgenes_neurons.pdf")

RUV1 <- RUVg(as.matrix(countData), controlGenes, k=1)
N1 <- RUV1$normalizedCounts

boxplot(log2(N1))
pca=princomp(cor(N1, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k1, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N1, colData = sampleData, design = ~ condition)
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



RUV2 <- RUVg(as.matrix(countData), controlGenes, k=2)
N2 <- RUV2$normalizedCounts

boxplot(log2(N2))
pca=princomp(cor(N2, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k2, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N2, colData = sampleData, design = ~ condition)
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



RUV3 <- RUVg(as.matrix(countData), controlGenes, k=3)
N3 <- RUV3$normalizedCounts

boxplot(log2(N3))
pca=princomp(cor(N3, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k3, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N3, colData = sampleData, design = ~ condition)
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



RUV4 <- RUVg(as.matrix(countData), controlGenes, k=4)
N4 <- RUV4$normalizedCounts

boxplot(log2(N4))
pca=princomp(cor(N4, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k4, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N4, colData = sampleData, design = ~ condition)
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


RUV5 <- RUVg(as.matrix(countData), controlGenes, k=5)
N5 <- RUV5$normalizedCounts

boxplot(log2(N5))
pca=princomp(cor(N5, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k5, neur proportion", pch=20, col=labels2colors(sampleData$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N5, colData = sampleData, design = ~ condition)
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



#removing samples <50
#aqui as amostras com propor��o neuronal menor do que 50% s�o removidas

countData50 <- counts [,c(1,2,5:7,9:15,17:23)]
sampleData50 <- sampleData [c(1,2,5:7,9:15,17:23), c(1:5)]
geneRPKM50 <- geneRPKM [,c(1,2,5:7,9:15,17:23)]

filtered <- which(rowSums(geneRPKM50 > 1) >= 9)
RPKM <- geneRPKM50[filtered,]
countData50 =subset (countData50, rownames(countData50) %in% rownames(RPKM))

#differential expression with Deseq (not considering batch in design formula)
Aqui� rodada a an�lise de diferen�a de express�o com essa lista de amostras >50% para gerar a listas de genes que n�o foram diferencialmente expressos (p>=0.8)

dds <- DESeqDataSetFromMatrix(countData = countData50, colData = sampleData50, design = ~ condition)

#set which group is the reference group
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
#you can check the number of outliers and lowcount genes using the function summary
summary(res)
res <- as.data.frame(res)
nDEGs_deseq50 <- subset (res, res$pvalue>=0.8)
nDEGs_deseq50 <- rownames(nDEGs_deseq50)

#control genes
HKgenes50 <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData50))

nDEGs_HKgenes50 <- intersect(nDEGs_deseq50,HKgenes50)
m=match(colnames(countData50), sampleData50$label)


#RUV

#before normalization
#Aqui os gr�ficos de PCA e de cluster hier�rquico s�o plotados com essa lista de amostras >50%

pdf("D:/backup_arquivos_estagio/ruv-seq_analysis/Results_RUVSeq_DROPBOX/BeforeNorm_plots_neurons50.pdf")
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


#nDEG_HKgenes
#Aqui o RUVseq � rodado considerando a lista de amostras >50% e a lista de nDEGs_HKgenes como a list de refer�ncia. Talvez aqui fa�a sentido mostrar os comandos do RUVseq somente para k=4 (como se tivessemos escolhido o k usando a lista com todas as amostras, n�o s� as > do que 50%)

controlGenes <- nDEGs_HKgenes50
pdf("D:/backup_arquivos_estagio/ruv-seq_analysis/Results_RUVSeq_DROPBOX/contrGenes_nDEGs_HKgenes_neurons50.pdf")


RUV4 <- RUVg(as.matrix(countData50), controlGenes, k=4)
N4 <- RUV4$normalizedCounts

boxplot(log2(N4))
pca=princomp(cor(N4, method="s"))
plot(pca$loadings[,1], pca$loadings[,2],xlab="PC1", ylab="PC2", main="PCA,RUV Normalisation k4, neur proportion", pch=20, col=labels2colors(sampleData50$Neur_prop_group[m]))
dds <- DESeqDataSetFromMatrix(countData = N4, colData = sampleData50, design = ~ condition)
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
