####Neuron Analysis - FINAL SCRIPT ###########

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);

## Load required packages ##########################################
library(edgeR)
library(RUVSeq)
library(EDASeq)
library(WGCNA)
library(calibrate)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(biomaRt)
# you might have to download additional packages depending of your currently R library

## Load input data #############################################
#load FPKM table to select genes in count table by FPKM>1
setwd("C:/Users/Karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files")
geneRPKM=read.delim ("allData_FPKM_renormalized_IV.txt", sep= "\t", header=TRUE)
rownames(geneRPKM)=geneRPKM[,1]
geneRPKM=geneRPKM[,c(41:63)]

#Load HKgenes lists
HKgenes_full <- read.csv ("HK_full_gene_list_ensembl_biomart.csv")


#load countdata and sample data
setwd("C:/Users/Karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/Neuronios")
counts=read.table("countdata_20M_neurons.txt", header=TRUE)
sampleInfo=read.delim("samplesheet_neurons.txt", sep="\t", header=TRUE)

###Removing samples with <50%######
#Our final analysis was made with samples only from with more than 50% of neuron proportion
countData50 <- counts [,c(1,2,5:7,9:15,17:23)]
sampleData50 <- sampleInfo [c(1,2,5:7,9:15,17:23),]
geneRPKM50 <- geneRPKM [,c(1,2,5:7,9:15,17:23)]


filtered <- which(rowSums(geneRPKM50 > 1) >= 9)
RPKM <- geneRPKM50[filtered,]
countData50 =subset (countData50, rownames(countData50) %in% rownames(RPKM))

#find the intersect between countData and HKgenes
HKgenes <- intersect(HKgenes_full$Ensembl_biomart, rownames(countData50))


####RUVseq normalization#####
#differential expression with Deseq
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

nDEGs_HKgenes <- intersect(nDEGs_deseq50,HKgenes)
controlGenes <- nDEGs_HKgenes

RUV4 <- RUVg(as.matrix(countData50), controlGenes, k=4)
N4 <- RUV4$normalizedCounts


write.table (N4, file="countData_RUV_Normalized_Neuron_nDEGs_HKgenes_p0.4_k4.txt")

###DESeq2 - Differential expression analysis#####

#create deseq object
dds <- DESeqDataSetFromMatrix(countData = N4, colData = sampleData50, design = ~ condition)

dds$condition <- relevel(dds$condition, ref="C")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
summary(res)
write.csv( as.data.frame(res), file="deseq_notCol_Neuron_nDEGs_HKgenes_p0.4_k4.csv" )


#to find out which is the outlier sample: download the following table
#the outlier genes are those with "NA" in both pvalue and padj columns
#you can check in the "cooks distance" table, which is the outlier sample
outliers <- assays (ddsdeseq)[["cooks"]]
write.csv( as.data.frame(outliers), file="outliers_notCol_Neuron_nDEGs_HKgenes_p0.4_k4.csv" )
#To extract the mean normalized counts
mnc <- assays(ddsdeseq)[["mu"]]
write.csv( as.data.frame(mnc), file="meannormalizedcounts_notCol_Neuron_nDEGs_HKgenes_p0.4_k4.csv" )

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
ncvsd2<- ncvsd
ncvsd2[,20]<-rownames (ncvsd2)
colnames(ncvsd2)[20]="Symbol"
resout <- as.data.frame (res)
resout [,7] <- rownames (resout)
colnames(resout)[7]="Symbol"
logtrans <- left_join(resout, ncvsd2, by = "Symbol")
write.csv(logtrans, file="logtransfcounts_notCol_Neuron_nDEGs_HKgenes_p0.4_k4.csv" )


#####WGCNA analysis######
#check for variance over the mean and possibliy remove some genes
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)

keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.0125)#I calculate the variance over the mean and tried different cuttoffs to see how many genes were still retained
ncvsd=ncvsd[keep,]


#transpose the table
dat=t(ncvsd)
infoData=rownames(ncvsd)



#Run WGCNA#
# run soft pick threshold
# Test a series of powers to which co-expression similarity is raised to calculate adjacency
powers=c(c(1:10), seq(from = 12, to = 20, by = 2))
sft = pickSoftThreshold(dat, powerVector = powers, verbose = 5, blockSize = 20000, networkType = "signed") # Signed gives an indication of positive vs negative correlations
#Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9 #what is that?
# Fit to scale-free topology
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2", type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels = powers, cex = cex1, col="red");
abline(h=c(0.8, 0.5), col="red")
# Mean connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab = "Soft threshold (power)",
     ylab = "Mean Connectivity", type = "n", main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col="red")
pdf(file="WGCNA_SoftThreshold_Neuron_nDEGs_HKgenes_p0.4_k4.pdf")
dev.off()

#construction of the network
# Chosen soft threshold in this case is 12; it sets the fit to Scale-free Topology > 0.8 while retaining as much connectivity as possible
net=blockwiseModules(dat, power=16, numericLabels=TRUE, networkType = "signed",
                     minModuleSize=150, mergeCutHeight=0.15, saveTOMs=FALSE, verbose=6,minKMEtoStay = 0.5,
                     nThreads=24, maxBlockSize=20000, checkMissingData=FALSE)


## Labelling the modules with a colour tag
modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N")
modules$Label=paste("M", modules$Label, sep="");
modules$Color=c("grey",labels2colors(modules$Label[-1]))
moduleLabel=paste("M",net$colors, sep="");
moduleColor=modules$Color[match(moduleLabel, modules$Label)]

####dendogram modules
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#####get kMe table####
## Calculating kMEs
KMEs<-signedKME(dat, net$MEs,outputColumnName = "M") # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme=data.frame(infoData [match(colnames(dat), infoData)], moduleColor,moduleLabel, KMEs)
colnames(kme)[1]="Symbol"

#Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedata=kme[,-kmeInfoCols];
pvalBH=kmedata; pvalBH[,]=NA
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedata[,j], nSamples=19), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme$newModule="NA"
for (j in c(1:nrow(kmedata)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedata)))
  m=which(kmedata[j,]==max(kmedata[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedata[j,m]>0.5)) kme$newModule[j]=as.character(colnames(kmedata)[m])
}

## Assign genes not associated to any module to M0
kme$newModule[which(kme$newModule%in%"NA")]="M0"

## Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme$newColor=kme$moduleColor[match(kme$newModule, kme$moduleLabel)]
kme$moduleLabel=kme$newModule; kme$moduleColor=kme$newColor
kme=kme[,-grep("newModule", colnames(kme))];kme=kme[,-grep("newColor", colnames(kme))]

#Saving kMEs
mod=modules$Label[-1]
kmeTable=kme[,kmeInfoCols];
for(j in c(1:length(mod)))
{
  kmeTable=cbind(kmeTable, kmedata[,match(mod[j],colnames(kmedata))]);colnames(kmeTable)[ncol(kmeTable)]=paste("kME", mod[j], sep="_")
  kmeTable=cbind(kmeTable, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTable)[ncol(kmeTable)]=paste("pvalBH", mod[j], sep="_")
}
# Renaming the rBEE rows, which are now NA
write.csv(kmeTable, "kME_Neuron_nDEGs_HKgenes_p0.4_k4.csv", row.names=FALSE)

#Saving Module Eigengenes
me<-data.frame(rownames(dat), net$MEs) # Sample names bound to module eigengenes
colnames(me)[-1]=gsub("ME", "M", colnames(me)[-1])
colnames(me)[1]="Sample"
write.csv(me, "ME_Neuron_nDEGs_HKgenes_p0.4_k4.csv", row.names=FALSE)

#Plotting Module Eigengene Values
colors=rep("turquoise", nrow(me) )
colors[grep("_P_", me$Sample)]="red"

pdf("moduleBarplots_Neuron_nDEGs_HKgenes_p0.4_k4.pdf", height=5, width=15)

mod=paste("M", c(0:(ncol(me)-2)), sep="")
for(m in mod)
{
  j=match(m, colnames(me))
  col=kme$moduleColor[match(m, kme$moduleLabel)]
  barplot(me[,j],  xlab="Samples", ylab="ME",col=colors, main=m, names=me[,1], cex.names=0.5, axisnames = FALSE)

}
dev.off()

# get biological variables correlations
# Define numbers of genes and samples
Samples = rownames(dat);
traitRows = match(Samples, sampleData50$label);
SampleInfo = sampleData50[traitRows, -1];
rownames(SampleInfo) = sampleData50[traitRows, 1];
nGenes = ncol(dat);
nSamples = nrow(dat);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dat, moduleColor)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, SampleInfo, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);



sizeGrWindow(10,6)
pdf("Module-TRait_relationship_Neuron_nDEGs_HKgenes_p0.4_k4.pdf")
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 2), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
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
               main = paste("Module-trait relationships"))
dev.off()
pvalue = moduleTraitPvalue [,3]
n = c(15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15)
bf = p.adjust(0.04, method = "bonferroni", n = 25)
bf

####overlapping of modules with external data####
#converting ensembl id to gene symbol
ensembl <- rownames (ncvsd)
geneR <-if (interactive()) {
  mart <- useMart("ensembl")
  datasets <- listDatasets(mart)
  mart <- useDataset("hsapiens_gene_ensembl",mart)
  getBM(attributes = c("hgnc_symbol", "ensembl_gene_id"),
        filters    = "ensembl_gene_id",
        values     = ensembl ,
        mart       = mart)
}


#comparison with brainLists and brainMarkers
colnames (geneR) <- c("Gene", "Symbol")
labelR <- left_join (geneR, kmeTable, by = "Symbol")
write.csv(labelR, "kME_Neuron_nDEGs_HKgenes_p0.4_k4.csv", row.names=FALSE)
labelR <- labelR [,3]
geneR <- geneR [,1]
data("BrainLists")
data ("BrainRegionMarkers")
userListEnrichment(
  geneR, labelR,
  fnIn = NULL,
  nameOut = "enrichment_brainLists_Neuron_nDEGs_HKgenes_p0.4_k4.csv",
  useBrainLists = TRUE, omitCategories = "grey",
  useBrainRegionMarkers = TRUE)

userListEnrichment(
  geneR, labelR,
  fnIn = NULL,
  nameOut = "enrichment_bloodLists_Neuron_nDEGs_HKgenes_p0.4_k4.csv",
  useBloodAtlases = TRUE, omitCategories = "grey")


mod_mari <- read.csv ("analises_50/modules_Mariani.csv")

userListEnrichment(
  geneR, labelR,
  fnIn = ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/analises_50/modules_Mariani.csv"),
  nameOut = "enrichment_mariani_Neuron_nDEGs_HKgenes_p0.4_k4.csv",
  omitCategories = "grey")


userListEnrichment(
  geneR, labelR,
  fnIn = ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/modules_Gupta.csv"),
  nameOut = "enrichment_Gupta_Neuron_nDEGs_HKgenes_p0.4_k4.csv",
  omitCategories = "grey")
