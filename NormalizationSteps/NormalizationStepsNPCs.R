###NormalizationStepsNPCs###

#display the current working directory
getwd();
#if necessary, change the path below to the directory where the data files are stored.
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
#you might have to download additional packages depending on your currently R library

#the following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#loading input data
geneFPKM=read.delim ("allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneFPKM)=geneFPKM[,1]
geneFPKM=geneFPKM[,-1]

counts=read.table("countdata_20M_NPC.txt", header=TRUE)
sampleInfo=read.delim("samplesheet_NPC.txt", sep="\t", header=TRUE)
HKgenes_full <- read.csv ("HK_full_gene_list_ensembl_biomart.csv")

#removing samples with >55%
#samples with neuronal cell proportion > 55% are removed
countData <- counts [,c(1:10,12:14,22:37)]
sampleData <- sampleInfo [c(1:10,12:14,22:37),]
sampleData <- sampleData[,c(1:2,4,9:10)]
geneFPKM2 <- geneFPKM [,c(1:10,12:14,22:37)]

filtered <- which(rowSums(geneFPKM2 > 1) >= 15)
FPKM <- geneFPKM2[filtered,]
countData=subset (countData, rownames(countData) %in% rownames(FPKM))

#find the intersect between countData and HKgenes
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData))


#first differential expression analysis with DESeq2
#we generated a gene list which is known to be no differentially expressed (pvalue>=0.7) so we can built up a housekeeping gene list for RUVseq analysis
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)

#set which group is the reference group
dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)

#you can check the number of outliers and lowcount genes using the function summary
summary(res)
res <- as.data.frame(res)
nDEGs_deseq <- subset (res, res$padj>=0.7)

#list of housekeeping genes used for the RUVseq normalization analysis
nDEGs_deseq <- rownames(nDEGs_deseq)
nDEGs_HKgenes <- intersect(nDEGs_deseq,HKgenes)
controlGenes <- nDEGs_HKgenes
m=match(colnames(countData), sampleData$label)

###RUVSeq###

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

#testing different K values which is related to the number of factors of unwanted variation
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
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition) 

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)


###WGCNA analysis###
#check for variance over the mean and possibliy remove some genes
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)

#we calculated the variance over the mean and tried different cuttoffs to see how many genes were still retained
keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.020)
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

pdf(file="BeforeNormalizationWGCNA_SoftThreshold_NPCs_nDEGs_HKgenes_p0_7.pdf")
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

#construction of the network
#chosen soft threshold in this case is 12; it sets the fit to Scale-free Topology > 0.8 while retaining as much connectivity as possible
net=blockwiseModules(dat, power=10, numericLabels=TRUE, networkType = "signed",
                     minModuleSize=50, mergeCutHeight=0.15, saveTOMs=FALSE, verbose=6,minKMEtoStay = 0.5,
                     nThreads=24, maxBlockSize=20000, checkMissingData=FALSE)

#labelling the modules with a colour tag
modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N")
modules$Label=paste("M", modules$Label, sep="");
modules$Color=c("grey",labels2colors(modules$Label[-1]))
moduleLabel=paste("M",net$colors, sep="");
moduleColor=modules$Color[match(moduleLabel, modules$Label)]

#get kMe table
#calculating kMEs
KMEs<-signedKME(dat, net$MEs,outputColumnName = "M") # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme=data.frame(infoData[match(colnames(dat), infoData)], moduleColor,moduleLabel, KMEs)
colnames(kme)[1]="Symbol"

#re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedata=kme[,-kmeInfoCols];
pvalBH=kmedata; pvalBH[,]=NA
for (j in c(1:ncol(pvalBH))){
  p=mt.rawp2adjp(corPvalueStudent(kmedata[,j], nSamples=29), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme$newModule="NA"
for (j in c(1:nrow(kmedata)) ){
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m=which(kmedata[j,]==max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedata[j,m]>0.5)) kme$newModule[j]=as.character(colnames(kmedata)[m])
}

#assign genes not associated to any module to M0
kme$newModule[which(kme$newModule%in%"NA")]="M0"

#replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)]
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor
kme=kme[,-grep("newModule", colnames(kme))];kme=kme[,-grep("newColor", colnames(kme))]

#saving kMEs
mod=modules$Label[-1]
kmeTable=kme[,kmeInfoCols];
for(j in c(1:length(mod))){
  kmeTable=cbind(kmeTable, kmedata[,match(mod[j],colnames(kmedata))]);colnames(kmeTable)[ncol(kmeTable)]=paste("kME", mod[j], sep="_")
  kmeTable=cbind(kmeTable, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTable)[ncol(kmeTable)]=paste("pvalBH", mod[j], sep="_")
}

#renaming the rBEE rows, which are now NA
write.csv(kmeTable, "BeforeNormalization_kME_NPC_nDEGs_HKgene_p0.7.csv", row.names=FALSE)

#saving module eigengenes
me<-data.frame(rownames(dat), net$MEs) # Sample names bound to module eigengenes
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me)[1]="Sample"
write.csv(me, "BeforeNormalization_ME_NPC_nDEGs_HKgene_p0.7.csv", row.names=FALSE)

#plotting module eigengene values
colors=rep("turquoise", nrow(me) )
colors[grep("_P_", me$Sample)]="red"

pdf("BeforeNormalization_moduleBarplots_NPC_nDEGs_HKgene_p0.7.pdf", height=5, width=15)

 mod=paste("M", c(0:(ncol(me)-2)), sep="")
 for(m in mod){
   j=match(m, colnames(me))
   col=kme$moduleColor[match(m, kme$moduleLabel)]
   barplot(me[,j],  xlab="Samples", ylab="ME",col=colors, main=m, names=me[,1], cex.names=0.5, axisnames = FALSE)
 }
dev.off()

#get biological variables correlations
#define numbers of genes and samples
Samples = rownames(dat);
traitRows = match(Samples, sampleData$label);
SampleInfo = sampleData[traitRows, -1];
rownames(SampleInfo) = sampleData[traitRows, 1];
nGenes = ncol(dat);
nSamples = nrow(dat);
#recalculate MEs with color labels
MEs0 = moduleEigengenes(dat, moduleColor)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, SampleInfo, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)

pdf("BeforeNormalization_Module-TRait_NPC_nDEG_HKgenes_p0.7.pdf")

 #will display correlations and their p-values
 textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
 dim(textMatrix) = dim(moduleTraitCor)
 par(mar = c(6, 8.5, 3, 3));

 #display the correlation values within a heatmap plot
 labeledHeatmap(Matrix = moduleTraitCor,
                xLabels = names(SampleInfo),
                yLabels = names(MEs),
                ySymbols = names(MEs),
                colorLabels = FALSE,
                colors = greenWhiteRed(50),
                textMatrix = textMatrix,
                setStdMargins = FALSE,
                cex.text = 0.5,
                zlim = c(-1,1),
                main = paste("BeforeNormalization_Module-trait relationships"))
 dev.off()

#once we have generated PCA and cluster hierarchical graphs considering K=1,..., K=4 in the RUVSeq and after having a look at the WGCNA results prior to normalization we were able to choose the best K value which would remove the observed batch effect related to neuronal cell proportions in our samples 

RUV2 <- RUVg(as.matrix(countData), controlGenes, k=2)
N2 <- RUV2$normalizedCounts
write.table (N2, file="countData_RUV_Normalized_NPC_nDEGs_HKgenes_p0.7_k2.txt")

#further steps after normalization are located at Final_script_NPC_analysis.R 