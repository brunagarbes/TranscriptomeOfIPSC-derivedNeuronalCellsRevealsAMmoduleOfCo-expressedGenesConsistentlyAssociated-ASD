###UserListEnrichment Neurons VS NPCs ###

#display the current working directory
getwd();
#if necessary, change the path below to the directory where the data files are stored.
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir);

#load required packages

library(stringr)
library(dplyr);
library(readr);
library(magrittr);
library(readxl);
library(WGCNA);
#you might have to download additional packages depending on your currently R library

#the following setting is important, do not omit.
options(stringsAsFactors = FALSE);

KMENeuron <- read_excel("kME_counts_Neurons_p0.4_k4_keep0.0125_power16_minsize150_cut0.15.xlsx")
KMENeuron <- select(KMENeuron, Symbol : moduleLabel)
colnames(KMENeuron) <- c("Ensembl","GeneName", "moduleColor", "moduleLabel")

KMENPC <- read_csv("kME_NPC_nDEGs_HKgene_p0.3_k2.csv")
KMENPC = select(KMENPC, Gene:moduleLabel)
colnames(KMENPC)= c("GeneName", "Ensembl", "moduleColor", "moduleLabel")


geneR=KMENeuron$Ensembl
labelR=KMENeuron$moduleColor
refTable=select(KMENPC, Ensembl, moduleColor)
write.csv(refTable, "refTable.csv", row.names = FALSE)

results = userListEnrichment(
  geneR, labelR, 
  fnIn = ("refTable.csv"),
  catNmIn = c("NPCRefList"),
  nameOut = "enrichmentBrainLists_NPC_fnIn_vs_Neuron.csv", 
  omitCategories = "grey",
  outputCorrectedPvalues = TRUE) 

results$sigOverlaps
head(results$pValue)
results$ovGenes

sigOverlap=results$sigOverlaps
sigOverlap=select(sigOverlap, InputCategories, UserDefinedCategories, CorrectedPvalues)
colnames(sigOverlap)= c("InputCategories_Neuron", "UserDefinedCategories_NPCrefList", "CorrectedPvalues")
UserDefinedCategories_NPC= str_replace(sigOverlap$UserDefinedCategories,"__NPCRefList", " ")
sigOverlap= cbind(sigOverlap, UserDefinedCategories_NPC)
sigOverlap=select(sigOverlap,InputCategories_Neuron, UserDefinedCategories_NPC, CorrectedPvalues)



geneR2=KMENPC$Ensembl
labelR2=KMENPC$moduleColor
refTable2=select(KMENeuron,Ensembl, moduleColor)
write.csv(refTable2, "reftable2.csv", row.names = FALSE)

results2 = userListEnrichment(
  geneR2, labelR2, 
  fnIn = ("reftable2.csv"),
  catNmIn = c("NeuronRefList"),
  nameOut = "enrichmentBrainLists_Neuron_fnIn_vs_NPC.csv", 
  omitCategories = "grey", 
  outputCorrectedPvalues = TRUE)

results2$sigOverlaps
head(results2$pValue)
results2$ovGenes

sigOverlap2=results2$sigOverlaps
sigOverlap2=select(sigOverlap2, InputCategories, UserDefinedCategories, CorrectedPvalues)
colnames(sigOverlap2)= c("InputCategories_NPC", "UserDefinedCategories_NeuronRefList", "CorrectedPvalues")
UserDefinedCategories_Neuron= str_replace(sigOverlap2$UserDefinedCategories,"__NeuronRefList", " ")
sigOverlap2= cbind(sigOverlap2, UserDefinedCategories_Neuron)
sigOverlap2=select(sigOverlap2,InputCategories_NPC, UserDefinedCategories_Neuron, CorrectedPvalues)

