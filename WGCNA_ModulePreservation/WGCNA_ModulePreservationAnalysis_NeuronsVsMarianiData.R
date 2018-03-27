
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
#workingDir = ".";
#setwd(workingDir);
#loading packages
library(WGCNA);
library(plyr);
library(dplyr);
library(readr);
library(magrittr);
library(data.table);
library(readxl);
#you might have to download additional packages depending on your currently R library

#the following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#our neuron data 
ruvNormDataNeuron <- read_delim("countData_RUV_Normalized_Neuron_nDEGs_HKgenes_p0.4_k4.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE);

KMENeuron <- read_excel("kME_counts_Neurons_p0.4_k4_keep0.0125_power16_minsize150_cut0.15.xlsx");

KMENeuron <- select(KMENeuron, Symbol : moduleLabel);
colnames(KMENeuron) <- c("Gene","GeneName", "moduleColor", "moduleLabel");
ourDataNeuron <- merge(KMENeuron, ruvNormDataNeuron, by = ("Gene"));


#reading all .txt files containing mariani's 11 days expression data 

larquivos<-list.files("data_comparison/GSE61476_marianiData/11_days/counts_11dias",full.names=TRUE)
arquivos <- lapply(larquivos, function(x) fread(file = x))

arquivos <- lapply(larquivos, function(x) {
  df <- read.delim(file = x,header = F)
  return(df)
})

dfarquivos = data.frame(arquivos)
mariani11Days = dfarquivos %>% select(V4, starts_with("V5")) 
colnames(mariani11Days) = c("geneID","GSM1505819sp","GSM1505825sp","GSM1505826sp","GSM1505827sp","GSM1505828sp","GSM1505829sp","GSM1505830sp","GSM1505839sp","GSM1505840sp","GSM1505842sp","GSM1505844sp","GSM1505846sp","GSM1505848sp","GSM1505852sp","GSM1505854sp","GSM1505855sp","GSM1505858sp","GSM1505860sp","GSM1505862sp")

moduleGenes = read_csv("data_comparison/GSE61476_marianiData/modules_Mariani.csv")
colnames(moduleGenes) = c("geneID", "module"); #we use such information so we can know which genes the original author actually used in her analysis
mariani11days =  merge(mariani11Days, moduleGenes, by = ("geneID"))


#reading all .txt files containing mariani's 31 days expression data 

larquivos31dias <-list.files("data_comparison/GSE61476_marianiData/31_days/counts_31dias",full.names=TRUE)
arquivos <- lapply(larquivos31dias, function(x) fread(file = x))

arquivos31dias <- lapply(larquivos31dias, function(x) {
  df <- read.delim(file = x,header = F)
  return(df)
})

dfarquivos31dias = data.frame(arquivos31dias)
mariani31Days = dfarquivos31dias %>% select(V4, starts_with("V5")) 
colnames(mariani31Days) = c("geneID","GSM1505820sp","GSM1505821sp","GSM1505831sp","GSM1505832sp","GSM1505833sp","GSM1505834sp","GSM1505835sp","GSM1505836sp","GSM1505841sp","GSM1505843sp","GSM1505845sp","GSM1505847sp","GSM1505849sp","GSM1505850sp","GSM1505851sp","GSM1505853sp","GSM1505856sp","GSM1505857sp","GSM1505859sp", "GSM1505861sp", "GSM1505863sp")
mariani31days =  merge(mariani31Days, moduleGenes, by = ("geneID"));
marianiMergedData = merge(mariani11days, mariani31days, by = ("geneID"));
marianiMergedData = marianiMergedData %>% select(everything(), -module.y, -module.x)


#Module preservation analysis Neurons VS Mariani's 11 AND 31 days data

#reading neuron data
dat0 = ourDataNeuron
names(dat0)
# this contains information on the genes
datSummaryNeuron=dat0[,c(1:4)] 
#the following data frame contains
#the gene expression data: columns are genes, rows are arrays (samples)
datExprNeuron <- t(dat0[,5:22])
no.samples <- dim(datExprNeuron)[[1]]
dim(datExprNeuron)
#set the columns names to probe names
colnames(datExprNeuron) = datSummaryNeuron$GeneName
colorsNeuron = dat0$moduleColor

#reading Mariani 11 AND 31 days Data
dataMerged = marianiMergedData
datExprMerged = t(dataMerged[,2:41]);
colnames(datExprMerged) = dataMerged$geneID
setLabels = c("Neuron", "MarianiMerged");
multiExpr = list(Neuron = list(data = datExprNeuron), marianiMerged = list(data = datExprMerged));
multiColor = list(Neuron = colorsNeuron);


#running Module Preservation Analysis Neurons vs Mariani's Merged data

 system.time( {
  mp = modulePreservation(multiExpr, multiColor,
                          referenceNetworks = 1,
                          nPermutations = 200,
                          randomSeed = 1,
                          quickCor = 0,
                          verbose = 3)
 } );

#save the results
save(mp, file = "modulePreservationNeuronMarianiMergedData.RData");

#alternatively, if the data has already been calculated before, load them from disk:
#load(file="modulePreservationNeuronMarianiMergedData.RData")

#analysis and display of module preservation results

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

#compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#plotting results

#module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
#leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
#text labels for points
text = modColors[plotMods];
#auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
#main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
#start the plot
sizeGrWindow(10, 5);

pdf(file ="NeuronsvsMarianiModulePreservationMergedData-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
 for (p in 1:2)
 {
   min = min(plotData[, p], na.rm = TRUE);
   max = max(plotData[, p], na.rm = TRUE);
   # Adjust ploting ranges appropriately
   if (p==2)
   {
     if (min > -max/10) min = -max/10
     ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
   } else
     ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
   plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
        main = mains[p],
        cex = 2.4,
        ylab = mains[p], xlab = "Module size", log = "x",
        ylim = ylim,
        xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
   labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
   # For Zsummary, add threshold lines
   if (p==2)
   {
     abline(h=0)
     abline(h=2, col = "blue", lty = 2)
     abline(h=10, col = "darkgreen", lty = 2)
   }
 }
dev.off();


#re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
#exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
#create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
#start the plot: open a suitably sized graphical window and set sectioning and margins.
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
#plot each Z statistic in a separate plot.
 for (s in 1:ncol(statsZ))
 {
   min = min(statsZ[plotMods, s], na.rm = TRUE);
   max = max(statsZ[plotMods, s], na.rm = TRUE);
   if (min > -max/5) min = -max/5
   plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
        main = colnames(statsZ)[s],
        cex = 1.7,
        ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
        xlim = c(20, 1000))
   labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
   abline(h=0)
   abline(h=2, col = "blue", lty = 2)
   abline(h=10, col = "darkgreen", lty = 2)
 }

data.frame(color = modColors[plotMods], label = labs)



#Module preservation analysis Neurons VS Mariani's 11 days data

#reading neuron data

dat0 = ourDataNeuron
names(dat0)
#this contains information on the genes
datSummaryNeuron=dat0[,c(1:4)] 
#the following data frame contains
#the gene expression data: columns are genes, rows are arrays (samples)
datExprNeuron <- t(dat0[,5:22])
no.samples <- dim(datExprNeuron)[[1]]
dim(datExprNeuron)
#set the columns names to probe names
colnames(datExprNeuron) = datSummaryNeuron$GeneName
colorsNeuron = dat0$moduleColor

#reading Mariani 11 days Data
dataEleven = mariani11days
datExprEleven = t(dataEleven[,2:20]);
colnames(datExprEleven) = dataEleven$geneID
setLabels = c("Neuron", "MarianiElevenDays");
multiExpr = list(Neuron = list(data = datExprNeuron), mariani11Days = list(data = datExprEleven));
multiColor = list(Neuron = colorsNeuron);

#running WGCNA module preservation analysis

 system.time( {
   mp = modulePreservation(multiExpr, multiColor,
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
 } );

#save the results
save(mp, file = "modulePreservationNeuronMarianiElevenDays.RData");

#alternatively, if the data has already been calculated before, load them from disk:
load(file="modulePreservationNeuronMarianiElevenDays.RData")


#analysis and display of module preservation results

ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

#compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#plotting results
#module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
#leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
#text labels for points
text = modColors[plotMods];
#auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
#main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
#start the plot
sizeGrWindow(10, 5);
pdf(file ="NeuronsvsMarianiModulePreservation11Days-Zsummary-medianRank.pdf", wi=10, h=5)
 par(mfrow = c(1,2))
 par(mar = c(4.5,4.5,2.5,1))
 for (p in 1:2)
 {
   min = min(plotData[, p], na.rm = TRUE);
   max = max(plotData[, p], na.rm = TRUE);
   # Adjust ploting ranges appropriately
   if (p==2)
   {
     if (min > -max/10) min = -max/10
     ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
   } else
     ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
   plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
        main = mains[p],
        cex = 2.4,
        ylab = mains[p], xlab = "Module size", log = "x",
        ylim = ylim,
        xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
   labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
   #for Zsummary, add threshold lines
   if (p==2)
   {
     abline(h=0)
     abline(h=2, col = "blue", lty = 2)
     abline(h=10, col = "darkgreen", lty = 2)
   }
 }

dev.off();

#re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
#exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
#create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
#plot each Z statistic in a separate plot.
 for (s in 1:ncol(statsZ))
 {
   min = min(statsZ[plotMods, s], na.rm = TRUE);
   max = max(statsZ[plotMods, s], na.rm = TRUE);
   if (min > -max/5) min = -max/5
   plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
        main = colnames(statsZ)[s],
        cex = 1.7,
        ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
        xlim = c(20, 1000))
   labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
   abline(h=0)
   abline(h=2, col = "blue", lty = 2)
   abline(h=10, col = "darkgreen", lty = 2)
 }

data.frame(color = modColors[plotMods], label = labs)


#Module preservation analysis Neurons VS Mariani's 31 days data

#reading neuron data
dat0 = ourDataNeuron
names(dat0)
#this contains information on the genes
datSummaryNeuron=dat0[,c(1:4)] 
#the following data frame contains
#the gene expression data: columns are genes, rows are arrays (samples)
datExprNeuron <- t(dat0[,5:22])
no.samples <- dim(datExprNeuron)[[1]]
dim(datExprNeuron)
#set the columns names to probe names
colnames(datExprNeuron) = datSummaryNeuron$GeneName
colorsNeuron = dat0$moduleColor

#reading Mariani 31 days Data
dataThirtyOne = mariani31days
#symbol = select(voiData, gene_symbol)
datExprThirtyOne = t(dataThirtyOne[,2:22]);
colnames(datExprThirtyOne) = dataThirtyOne$geneID
setLabels = c("Neuron", "MarianiThirtyOneDays");
multiExpr = list(Neuron = list(data = datExprNeuron), mariani31Days = list(data = datExprThirtyOne));
multiColor = list(Neuron = colorsNeuron);


#running WGCNA module preservation analysis

 system.time( {
   mp = modulePreservation(multiExpr, multiColor,
                           referenceNetworks = 1,
                           nPermutations = 200,
                           randomSeed = 1,
                           quickCor = 0,
                           verbose = 3)
 } );

#save the results
save(mp, file = "modulePreservationNeuronMarianiThirtyOneDays.RData");

#alternatively, if the data has already been calculated before, load them from disk:
#load(file="modulePreservationNeuronMarianiThirtyOneDays.RData")


#analysis and display of module preservation results
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);

#compare preservation to quality:
ZscoreTable=print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                         signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

#plotting results
#module labels and module sizes are also contained in the results
modColors = rownames(mp$preservation$observed[[ref]][[test]])
moduleSizes = mp$preservation$Z[[ref]][[test]][, 1];
#leave grey and gold modules out
plotMods = !(modColors %in% c("grey", "gold"));
#text labels for points
text = modColors[plotMods];
#auxiliary convenience variable
plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
#main titles for the plot
mains = c("Preservation Median rank", "Preservation Zsummary");
#start the plot
sizeGrWindow(10, 5);

pdf(file ="NeuronsvsMariani31daysModulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
 par(mfrow = c(1,2))
 par(mar = c(4.5,4.5,2.5,1))
 for (p in 1:2)
 {
   min = min(plotData[, p], na.rm = TRUE);
   max = max(plotData[, p], na.rm = TRUE);
   # Adjust ploting ranges appropriately
   if (p==2)
   {
     if (min > -max/10) min = -max/10
     ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
   } else
     ylim = c(max + 0.1 * (max-min), min - 0.1 * (max-min))
   plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
        main = mains[p],
        cex = 2.4,
        ylab = mains[p], xlab = "Module size", log = "x",
        ylim = ylim,
        xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4)
   labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08);
   # For Zsummary, add threshold lines
   if (p==2)
   {
     abline(h=0)
     abline(h=2, col = "blue", lty = 2)
     abline(h=10, col = "darkgreen", lty = 2)
   }
 }

dev.off();

#re-initialize module color labels and sizes
modColors = rownames(statsZ)
moduleSizes = mp$quality$Z[[ref]][[test]][, 1];
#exclude improper modules
plotMods = !(modColors %in% c("grey", "gold"));
#create numeric labels for each module
labs = match(modColors[plotMods], standardColors(50));
sizeGrWindow(12, 9);
par(mfrow = c(3,5))
par(mar = c(3,3,2,1))
par(mgp = c(1.6, 0.4, 0));
#plot each Z statistic in a separate plot.
 for (s in 1:ncol(statsZ))
 {
   min = min(statsZ[plotMods, s], na.rm = TRUE);
   max = max(statsZ[plotMods, s], na.rm = TRUE);
   if (min > -max/5) min = -max/5
   plot(moduleSizes[plotMods], statsZ[plotMods, s], col = 1, bg = modColors[plotMods], pch = 21,
        main = colnames(statsZ)[s],
        cex = 1.7,
        ylab = colnames(statsZ)[s], xlab = "Module size", log = "x",
        ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min)),
        xlim = c(20, 1000))
   labelPoints(moduleSizes[plotMods], statsZ[plotMods, s], labs, cex = 0.7, offs = 0.04);
   abline(h=0)
   abline(h=2, col = "blue", lty = 2)
   abline(h=10, col = "darkgreen", lty = 2)
 }

data.frame(color = modColors[plotMods], label = labs)
