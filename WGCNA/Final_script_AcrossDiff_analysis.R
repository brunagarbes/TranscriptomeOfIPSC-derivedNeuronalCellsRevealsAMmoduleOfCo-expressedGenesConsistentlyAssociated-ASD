###Differentiation Analysis - FINAL SCRIPT######
#Analysis of differential expression between NPC and Neuron, and which genes are differentially
#regulated between patients and controls across differentiation
##Load libraries####
library(multtest)
library(gplots)
library (DESeq2)
library(WGCNA)
library("RColorBrewer")
library("pheatmap")
library(RUVSeq)
library (dplyr)
library(biomaRt)
options(stringsAsFactors = FALSE)

####Loading Data#####
setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/NPC")
NPCData <- read.table ("countData_RUV_Normalized_NPC_nDEGs_HKgenes_p0.7_k2.txt", header=TRUE)
colnames (NPCData) [1] = "Gene"


setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/Neuronios")
NeuronData <- read.table ("countData_RUV_Normalized_Neuron_nDEGs_HKgenes_p0.4_k4.txt", header = TRUE)
colnames (NeuronData) [1]= "Gene"

countDataA = full_join (NPCData,NeuronData, by="Gene")
write.table (countDataA, file="countData_normalized_provisorio.txt")

setwd("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files/Resultados_finais/NPCxNeurons")
sampleData <- read.table ("sample_sheet_all_cells.txt", header = TRUE)
sampleData <- sampleData [c(1:10,12:14,22:37,40, 41, 44, 45, 46, 48:54, 56:62),]

###DESeq2 - Differential expression analysis#####

#create deseq object
dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ cell)

dds$cell <- relevel(dds$cell, ref="NPC")
ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq)
summary(res)
write.csv( as.data.frame(res), file="deseq_notCol_NPC_versus_neurons.csv" )


#to find out which is the outlier sample: download the following table
#the outlier genes are those with "NA" in both pvalue and padj columns
#you can check in the "cooks distance" table, which is the outlier sample
outliers <- assays (ddsdeseq)[["cooks"]]
write.csv( as.data.frame(outliers), file="outliers_notCol_NPC_Neurons.csv" )
#To extract the mean normalized counts
mnc <- assays(ddsdeseq)[["mu"]]
write.csv( as.data.frame(mnc), file="meannormalizedcounts_notCol_NPC_Neurons.csv" )

#to extract the log-transformed counts
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
ncvsd2<- ncvsd
ncvsd2[,49]<-rownames (ncvsd2)
colnames(ncvsd2)[49]="Symbol"
resout <- as.data.frame (res)
resout [,7] <- rownames (resout)
colnames(resout)[7]="Symbol"
logtrans <- left_join(resout, ncvsd2, by = "Symbol")
write.csv(logtrans, file="logtransfcounts_notCol_NPC_Neurons.csv" )


###interaction formula#######

dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition + cell + condition:cell)

dds$cell <- relevel(dds$cell, ref="NPC")
dds$condition <- relevel (dds$condition, ref="C")

ddsdeseq <- DESeq(dds)
res <- results(ddsdeseq, name="conditionP.cellNeur")
summary(res)
write.csv( as.data.frame(res), file="deseq_notCol_NPC_Neurons_C_diff_P.csv" )

resControls<- results (ddsdeseq, contrast=c("cell","Neur","NPC"))
summary(resControls)
write.csv( as.data.frame(resControls), file="deseq_notCol_NPC_Neurons_Ctrl.csv" )

resPat <- results (ddsdeseq, list (c("cell_Neur_vs_NPC", "conditionP.cellNeur")))
summary(resPat)
write.csv( as.data.frame(resPat), file="deseq_notCol_NPC_Neurons_Pat.csv" )


#####WGCNA#####
#Preparing the data

dds <- DESeqDataSetFromMatrix(countData = countData, colData = sampleData, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="C")
###log-transform the data
vsd<- varianceStabilizingTransformation(dds)
ncvsd <- assay (vsd)
ncvsd <- as.data.frame(ncvsd)
variance <- apply(ncvsd,1,sd)/apply(ncvsd,1,mean)
hist(variance)
keep=which(apply(ncvsd,1,sd)/apply(ncvsd,1,mean)>= 0.025)#I calculate the variance over the mean and tried different cuttoffs to see how many genes were still retained
ncvsd=ncvsd[keep,]


# We work with two sets:
NPCData <- ncvsd [,c(1:29)]
NeuronData <- ncvsd [,c(30:48)]


nSets = 2;
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("NPC", "Neuron")
shortLabels = c("NPC", "Neuron")
# Form multi-set expression data: 
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = as.data.frame(t(NPCData)))
names(multiExpr[[1]]$data) = rownames (NPCData)
rownames(multiExpr[[1]]$data) = names(NPCData)
multiExpr[[2]] = list(data = as.data.frame(t(NeuronData)))
names(multiExpr[[2]]$data) = rownames(NeuronData)
rownames(multiExpr[[2]]$data) = names(NeuronData)
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)


#cluster tree
sampleTrees = list()
for (set in 1:nSets)
{
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
pdf(file = "SampleClustering_NPC_Neuron.pdf", width = 12, height = 12);
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7)
dev.off()

#sample info

SampleInfo = vector(mode="list", length = nSets)
for (set in 1:nSets)
{
  setSamples = rownames(multiExpr[[set]]$data)
  SampleRows = match(setSamples, sampleData$label)
  SampleInfo[[set]] = list(data = sampleData [SampleRows, -1])
  rownames(SampleInfo[[set]]$data) = sampleData[SampleRows, 1]
}
collectGarbage()
# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples
infoData=rownames(NPCData)
save(multiExpr, SampleInfo, nGenes, nSamples, setLabels, shortLabels, exprSize,
     file = "Consensus-dataInput.RData");

###Construction of the network####

#load data
setwd ("C:/Users/karin/Dropbox/Arquivos_genomica_autistas/RNAseq/RNAseq_files")
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData");


# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 5)[[2]], blockSize = 10000, networkType = "signed");
collectGarbage();

# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
pdf(file = "scaleFreeAnalysis_WGCNA_consensus_NPC_Neurons.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()


##constructition of network


net = blockwiseConsensusModules(
  multiExpr, power = 12, minModuleSize = 50, networkType = "signed", deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  saveTOMs = FALSE, verbose = 5, maxBlockSize=20000)

#store the module eingengene in the variable consMEs
consMEs = net$multiMEs;
# Convert the numeric labels to color labels
#moduleLabels = net$colors;
#moduleColors = labels2colors(moduleLabels)


modules=as.data.frame(table(net$colors)); colnames(modules)=c("Label", "N") 
modules$Label=paste("kMEta.ME", modules$Label, sep=""); 
modules$Color=c("grey",labels2colors(modules$Label[-1])) 
moduleLabel=paste("kMEta.ME",net$colors, sep=""); 
moduleColor=modules$Color[match(moduleLabel, modules$Label)] 

#check the modules generated
consTree = net$dendrograms[[1]];
sizeGrWindow(8,6);
pdf(file = "ConsensusDendrogram-auto_all_cells_50_indepFPKMfilter_p0.5_k5.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColor,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()

#### Calculating kMEs#####
set = 1
KMEs_NPC<-signedKME(multiExpr[[set]], consMEs[[1]]) # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme_NPC=data.frame(infoData [match(colnames(multiExpr[[set]]$data), infoData)], moduleColor,moduleLabel, KMEs_NPC)
colnames(kme_NPC)[1]="Symbol"
#Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedataNPC=kme_NPC[,-kmeInfoCols]; 
pvalBH=kmedataNPC; pvalBH[,]=NA 
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedataNPC[,j], nSamples=29), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme_NPC$newModule="NA" 
for (j in c(1:nrow(kmedataNPC)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedataNPC)))
  m=which(kmedataNPC[j,]==max(kmedataNPC[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedataNPC[j,m]>0.5)) kme_NPC$newModule[j]=as.character(colnames(kmedataNPC)[m])
}

## Assign genes not associated to any module to M0
kme_NPC$newModule[which(kme_NPC$newModule%in%"NA")]="M0" 

## Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme_NPC$newColor=kme_NPC$moduleColor[match(kme_NPC$newModule, kme_NPC$moduleLabel)] 
kme_NPC$moduleLabel=kme_NPC$newModule; kme_NPC$moduleColor=kme_NPC$newColor 
kme_NPC=kme_NPC[,-grep("newModule", colnames(kme_NPC))];kme_NPC=kme_NPC[,-grep("newColor", colnames(kme_NPC))] 

#Saving kMEs 
mod=modules$Label[-1] 
kmeTableNPC=kme_NPC[,kmeInfoCols]; 
for(j in c(1:length(mod)))
{
  kmeTableNPC=cbind(kmeTableNPC, kmedataNPC[,match(mod[j],colnames(kmedataNPC))]);colnames(kmeTableNPC)[ncol(kmeTableNPC)]=paste("kME", mod[j], sep="_")                                                                             
  kmeTableNPC=cbind(kmeTableNPC, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTableNPC)[ncol(kmeTableNPC)]=paste("pvalBH", mod[j], sep="_")
}
# Renaming the rBEE rows, which are now NA    
write.csv(kmeTableNPC, "kME_consensus_NPC_var0.035_mod50_cut0.15.csv", row.names=FALSE)

set = 2
KMEs_Neurons<-signedKME(multiExpr[[set]], consMEs[[2]]) # signedKME calculates eigengene-based connectivity i.e. module membership. Also, outputColumnName gives us a prefix. Also, the corFnc defaults to Pearson
kme_Neurons=data.frame(infoData [match(colnames(multiExpr[[set]]$data), infoData)], moduleColor,moduleLabel, KMEs_Neurons)
colnames(kme_Neurons)[1]="Symbol"
#Re-assign genes to modules by kME. Any not passing the criteria of correlation pvalue <0.05, and kME>0.5, will be moved to a junk module
kmeInfoCols=c(1:3)
kmedataNeurons=kme_Neurons[,-kmeInfoCols]; 
pvalBH=kmedataNeurons; pvalBH[,]=NA 
for (j in c(1:ncol(pvalBH)))
{
  p=mt.rawp2adjp(corPvalueStudent(kmedataNeurons[,j], nSamples=19), proc="BH")
  pvalBH[,j]=p$adjp[order(p$index),2]
}

kme_Neurons$newModule="NA" 
for (j in c(1:nrow(kmedataNeurons)) )
{
  if (j==1) print("Working on genes 1:10000"); if(j==10000) print(paste("Working on genes 10000:", nrow(kmedataNeurons)))
  m=which(kmedataNeurons[j,]==max(kmedataNeurons[j,]))
  if ((pvalBH[j,m]<0.05)&(kmedataNeurons[j,m]>0.5)) kme_Neurons$newModule[j]=as.character(colnames(kmedataNeurons)[m])
}

## Assign genes not associated to any module to M0
kme_Neurons$newModule[which(kme_Neurons$newModule%in%"NA")]="M0" 

## Replacing the old moduleLabel and moduleColor columns with the updated newModule and newColor columns respectively
kme_Neurons$newColor=kme_Neurons$moduleColor[match(kme_Neurons$newModule, kme_Neurons$moduleLabel)] 
kme_Neurons$moduleLabel=kme_Neurons$newModule; kme_Neurons$moduleColor=kme_Neurons$newColor 
kme_Neurons=kme_Neurons[,-grep("newModule", colnames(kme_Neurons))];kme_Neurons=kme_Neurons[,-grep("newColor", colnames(kme_Neurons))] 

#Saving kMEs 
mod=modules$Label[-1] 
kmeTableNeurons=kme_Neurons[,kmeInfoCols]; 
for(j in c(1:length(mod)))
{
  kmeTableNeurons=cbind(kmeTableNeurons, kmedataNeurons[,match(mod[j],colnames(kmedataNeurons))]);colnames(kmeTableNeurons)[ncol(kmeTableNeurons)]=paste("kME", mod[j], sep="_")                                                                             
  kmeTableNeurons=cbind(kmeTableNeurons, pvalBH[,match(mod[j],colnames(pvalBH))]);colnames(kmeTableNeurons)[ncol(kmeTableNeurons)]=paste("pvalBH", mod[j], sep="_")
}
# Renaming the rBEE rows, which are now NA    
write.csv(kmeTableNeurons, "kME_consensus_Neurons_var0.035_mod50_cut0.15.csv", row.names=FALSE)


#Saving Module Eigengenes 
set = 1
me_NPC<-data.frame(rownames(multiExpr[[set]]$data), consMEs[[1]]) # Sample names bound to module eigengenes
colnames(me_NPC)[-1]=gsub("ME", "M", colnames(me_NPC)[-1])
colnames(me_NPC)[1]="Sample"
write.csv(me_NPC, "analises_50/ME_all_cells_consensus_NPC_50_indFPKMfilter_p0.5_k5_var0.035_mod50.csv", row.names=FALSE)

set = 2
me_neurons<-data.frame(rownames(multiExpr[[set]]$data), consMEs[[2]]) # Sample names bound to module eigengenes
colnames(me_neurons)[-1]=gsub("ME", "M", colnames(me_neurons)[-1])
colnames(me_neurons)[1]="Sample"
write.csv(me_neurons, "analises_50/ME_all_cells_consensus_neurons_50_indFPKMfilter_p0.5_k5_var0.035_mod50.csv", row.names=FALSE)


#Plotting Module Eigengene Values
set = 1
colors=rep("turquoise", nrow(me_Neurons) )
colors[grep("_NPC_P_", me_Neurons$Sample)]="red"

pdf("moduleBarplots_all_cells_consensus.pdf", height=5, width=15)

mod=paste("M", c(0:(ncol(me_Neurons)-2)), sep="")
for(m in mod) 
{ 
  j=match(m, colnames(me_Neurons))
  col=kme_NPC$moduleColors[match(m, kme_NPC$moduleLabels)]
  barplot(me_NPC[,j],  xlab="Samples", ylab="ME",col=colors, main=m, names=me_NPC[,1], cex.names=0.5, axisnames = FALSE)
  
}
dev.off()



#correlating the modules with sample traits
moduleTraitCor = list();
moduleTraitPvalue = list();
# Calculate the correlations
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, SampleInfo[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");
# Open a suitably sized window (the user should change the window size if necessary)
sizeGrWindow(10,7)
pdf(file = "analises_50/ModuleTraitRelationships_NPC_50_indFPKMfilter_p0.5_k5_var0.035_mod50.pdf", wi = 10, he = 7);
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(SampleInfo[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)
pdf(file = "analises_50/ModuleTraitRelationships_Neurons_50indFPKMfilter_p0.5_k5_var0.035_mod50.pdf", wi = 10, he = 7);
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(SampleInfo[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off()

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
labelR <- left_join (geneR, kmeTableNeurons, by = "Symbol")
labelR <- labelR [,3]
geneR <- geneR [,1]
data("BrainLists")
data ("BrainRegionMarkers")
userListEnrichment(
  geneR, labelR, 
  fnIn = NULL,
  nameOut = "enrichment_brainLists_consensus_analysis.csv", 
  useBrainLists = TRUE, omitCategories = "grey", 
  useBrainRegionMarkers = TRUE) 

mod_mari <- read.csv ("analises_50/modules_Mariani.csv")

userListEnrichment(
  geneR, labelR, 
  fnIn = ("analises_50/modules_Mariani.csv"),
  nameOut = "enrichment_brainLists_consensus_analysis_mariani.csv",
  omitCategories = "grey") 


