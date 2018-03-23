# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);
#loading packages
library(WGCNA);
library(plyr);
library(dplyr);
library(readr);
library(magrittr);
library(data.table);
library(readxl);
#the following setting is important, do not omit.
options(stringsAsFactors = FALSE);

#our neuron data 
ruvNormDataNeuron <- read_delim("countData_RUV_Normalized_Neuron_nDEGs_HKgenes_p0.4_k4.txt", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)

KMENeuron <- read_excel("kME_counts_Neurons_p0.4_k4_keep0.0125_power16_minsize150_cut0.15.xlsx")

KMENeuron <- select(KMENeuron, Symbol : moduleLabel)
colnames(KMENeuron) <- c("Gene","GeneName", "moduleColor", "moduleLabel")
ourDataNeuron <- merge(KMENeuron, ruvNormDataNeuron, by = ("Gene"))

#our NPC data

ruvNormDataNPC <- read_delim("data_comparison/NPC/countData_RUV_Normalized_NPC_nDEGs_HKgenes_p0.3_k2.txt", 
                             " ", escape_double = FALSE, trim_ws = TRUE)

KMENPC <- read_csv("data_comparison/NPC/kME_NPC_nDEGs_HKgene_p0.3_k2.csv")

KMENPC <- select(KMENPC, Symbol : moduleLabel)
colnames(KMENPC) <- c("Ensembl", "moduleColor", "moduleLabel")
ourDataNPC <- merge(KMENPC, ruvNormDataNPC, by = ("Ensembl"))


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
colnames(datExprNeuron) = datSummaryNeuron$Gene
#this module assignment was obtained by Ghazalpour et al
colorsNeuron = dat0$moduleColor

#reading NPC data

dataNPC = ourDataNPC
datExprNPC = t(dataNPC[,4:32]);
colnames(datExprNPC) = dataNPC$Ensembl
setLabels = c("Neuron", "NPC");
multiExpr = list(Neuron = list(data = datExprNeuron), NPC = list(data = datExprNPC));
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
save(mp, file = "modulePreservationNeuronNPC.RData");
#alternatively, if the data has already been calculated before, load them from disk:
#load(file="modulePreservationNeuronNPC.RData")

#analysis and display of module preservation results
ref = 1
test = 2
statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1]);
#compare preservation to quality:
print( cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
             signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )

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
mains = c("Preservation Median rank", "Preservation Zsummary - Neurons vs NPCs");
#start the plot
sizeGrWindow(10, 5);

pdf(file ="NeuronVsNPC-modulePreservation-Zsummary-medianRank.pdf", wi=10, h=5)
par(mfrow = c(1,2))
par(mar = c(4.5,4.5,2.5,1))
 for (p in 1:2)
 {
   min = min(plotData[, p], na.rm = TRUE);
   max = max(plotData[, p], na.rm = TRUE);
   #adjust ploting ranges appropriately
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