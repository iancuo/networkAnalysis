library(foreach)
library(doMC)
registerDoMC()
library(multtest)
library(WGCNA)
library("org.Mm.eg.db")
library(biomaRt)
library(GOstats)
library("org.Mm.eg.db")
library("edgeR")
library(vegan)
library(ncf)

library(lawstat)

#source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
 #biocLite("GOstats")

getDoParWorkers()
options(cores=27)
getDoParWorkers()

#enableWGCNAThreads(nThreads = 27)

setwd("/home/dan/workDir/networkAnalysis")

source("/home/dan/workDir/networkAnalysis/functionDefinitions.R")

load("selectedData_CEA.RData")

geneNames=rownames(selectedGeneCounts)

# divide the data in different groups
HSCC_H=selectedGeneCounts[,samplesHigh ]
HSCC_L=selectedGeneCounts[,samplesLow ]

results_highCPMgenes=read.csv("resultsCoexpr_CEA/resultsDEDV_highCPM.csv")

setdiff(rownames(selectedGeneCounts), results_highCPMgenes[,1])

########################################################################################################################################
adjConsensus=adjacency(t(selectedGeneCounts), power=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold.fromSimilarity(adjConsensus, powerVector = powers, verbose = 5, moreNetworkConcepts=T)

plotNetConstruction(sft)
quartz.save("figuresCoexpr/netConstructionCoexpr.tif", type="tif", bg="white", dpi=300)
quartz.save("figuresCoexpr/netConstructionCoexpr.jpg", type="jpg", bg="white")

softPowerCoexpr=6
adjCoexpr=adjConsensus^softPowerCoexpr
adjCoexpr[is.na(adjCoexpr)]=0

connCoexpr=rowSums(adjCoexpr)
hist(connCoexpr, 100)


hierADJCoexpr = hclust(as.dist(1-adjCoexpr),method="average");

# code below might need modifications to make the number of modules between ~ 15 and 50, and to make number of grey genes less than ~ 2-3000
hybridCoexpr=cutreeHybrid(dendro = hierADJCoexpr, distM=1-adjCoexpr, cutHeight = 0.9995, minClusterSize = 100, deepSplit = 4, maxCoreScatter = NULL, minGap = NULL, maxAbsCoreScatter = NULL, minAbsGap = NULL, pamStage = TRUE, pamRespectsDendro = F, useMedoids = FALSE,  respectSmallClusters = TRUE, verbose = 2, indent = 0)

colorsCoexpr = labels2colors(hybridCoexpr$labels)
names(colorsCoexpr)=geneNames
table(colorsCoexpr)
length(table(colorsCoexpr))
modulesCoexpr=names(table(colorsCoexpr))
sum(colorsCoexpr=="grey")


fileConnSummary<-file("resultsCoexpr_CEA/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Number of genes selected for network construction: ", length(geneNames), sep=''), fileConnSummary)

writeLines(paste("Number modules: ", length(table(colorsCoexpr)), sep=''), fileConnSummary)
writeLines(paste("Number grey genes:  ",  sum(colorsCoexpr=="grey"), sep=''), fileConnSummary)

close(fileConnSummary)

adj_HSCC_H=adjacency(t(HSCC_H), power=softPowerCoexpr)
adj_HSCC_L=adjacency(t(HSCC_L), power=softPowerCoexpr)

names(colorsCoexpr)=geneNames

save(adjCoexpr, colorsCoexpr,modulesCoexpr, adj_HSCC_H,adj_HSCC_L, file="data/adjCoexprModules.RData")
#load("data/adjCoexprModules.RData")
########################################################################################################
coexprConn=intramodularConnectivity(adjMat = adjCoexpr, colors = colorsCoexpr, scaleByMax = T)
coexprConnHigh=intramodularConnectivity(adjMat = adj_HSCC_H, colors = colorsCoexpr, scaleByMax = F)
coexprConnLow=intramodularConnectivity(adjMat = adj_HSCC_L, colors = colorsCoexpr, scaleByMax = F)

##############################################################
########################################################################################################

neuronsList=read.csv("data/CahoyNeurons.csv", header=TRUE)
neuronsSymbols= neuronsList[,"Gene.Name"]

astrosList=read.csv("data/CahoyAstros.csv", header=TRUE)
astrosSymbols= astrosList[,"Gene.Name"]

oligosList=read.csv("data/CahoyOligos.csv", header=TRUE)
oligosSymbols= oligosList[,"Gene.Name"]

moduleEnrichmentNeurons = moduleEnrichment (colorsCoexpr, neuronsSymbols)
moduleEnrichmentAstros = moduleEnrichment (colorsCoexpr, astrosSymbols)
moduleEnrichmentOligos = moduleEnrichment (colorsCoexpr, oligosSymbols)

cellTypeEnrichment=round(cbind(moduleEnrichmentNeurons,moduleEnrichmentAstros, moduleEnrichmentOligos),4)
colnames(cellTypeEnrichment)=c("Neurons", "Astros", "Oligos")
rownames(cellTypeEnrichment)=modulesCoexpr
write.csv(cellTypeEnrichment, file="resultsCoexpr_CEA/cellTypeEnrich.csv", append=T)

modulesEnrichedNeurons=modulesCoexpr[cellTypeEnrichment[,"Neurons"] < (0.05/length(modulesCoexpr))]
modulesEnrichedAstros=modulesCoexpr[cellTypeEnrichment[,"Astros"] < (0.05/length(modulesCoexpr))]
modulesEnrichedOligos=modulesCoexpr[cellTypeEnrichment[,"Oligos"] < (0.05/length(modulesCoexpr))]

modulesNeuros=""
for(i in 1:length(modulesEnrichedNeurons)){
  modulesNeuros=paste(modulesNeuros, modulesEnrichedNeurons[i])
}

modulesAstros=""
for(i in 1:length(modulesEnrichedAstros)){
  modulesAstros=paste(modulesAstros, modulesEnrichedAstros[i])
}

modulesOligos=""
for(i in 1:length(modulesEnrichedOligos)){
  modulesOligos=paste(modulesOligos, modulesEnrichedOligos[i])
}

fileConnSummary<-file("resultsCoexpr_CEA/SummaryResultsCoexpr.txt",  open="at")
writeLines(paste("Modules enriched in neuronal markers: ", modulesNeuros, sep=""), fileConnSummary)
writeLines(paste("Modules enriched in astrocyte markers: ", modulesAstros, sep=" "), fileConnSummary)
writeLines(paste("Modules enriched in oligodendrocyte markers: ", modulesOligos, sep=" "), fileConnSummary)

close(fileConnSummary)



########################################################################################################

# consensus module preservation in individual networks
# load("data/adjCoexprModules.RData")
# 
# multiData =vector("list",3)
# 
# multiData[[1]] =list(data= adjCoexpr)
# multiData[[2]] =list(data= adj_HSCC_H)
# multiData[[3]] =list(data= adj_HSCC_L)
# 
# names(multiData)=c("Consensus","HSCC_H", "HSCC_L")
# checkSets(multiData, checkStructure = FALSE, useSets = NULL)
# 
# multiColor =vector("list",3)
# 
# multiColor[[1]] =as.vector(colorsCoexpr)
# multiColor[[2]] =as.vector(colorsCoexpr)
# multiColor[[3]] =as.vector(colorsCoexpr)
# 
# names(multiColor)=c("Consensus","HSCC_H", "HSCC_L")
# 
# modulePreservIndVsConsensus=modulePreservation(
#   multiData,
#   multiColor,
#   dataIsExpr = F,
#   networkType = "unsigned", 
#   corFnc = "cor",
#   corOptions = "use = 'p'",
#   referenceNetworks = 1, 
#   nPermutations = 200, 
#   includekMEallInSummary = FALSE,
#   restrictSummaryForGeneralNetworks = FALSE,
#   calculateQvalue = FALSE,
#   randomSeed = 12345, 
#   maxGoldModuleSize = 1000, 
#   maxModuleSize = 1000, 
#   quickCor = 1, 
#   ccTupletSize = 2, 
#   calculateCor.kIMall = TRUE,
#   useInterpolation = FALSE, 
#   checkData = F, 
#   greyName = "grey", 
#   savePermutedStatistics = FALSE, 
#   loadPermutedStatistics = FALSE, 
#   permutedStatisticsFile = if (useInterpolation) "permutedStats-intrModules.RData" 
#   else "permutedStats-actualModules.RData", 
#   plotInterpolation = FALSE, 
#   interpolationPlotFile = "modulePreservationInterpolationPlots.pdf", 
#   discardInvalidOutput = TRUE,
#   verbose = 3, indent = 0)
# 
# 
# save(modulePreservIndVsConsensus, file="data/modulePreservation.RData")
# load("data/modulePreservation.RData")
# 
# names(modulePreservIndVsConsensus)
# names(modulePreservIndVsConsensus$preservation)
# names(modulePreservIndVsConsensus$preservation$Z)
# names(modulePreservIndVsConsensus$preservation$Z$ref.Consensus)
# names(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H)
# 
# 
# preservSummary=cbind(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H$moduleSize,
#                      modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H$Zsummary.pres,
#                      modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_L$Zsummary.pres)
# 
# rownames(preservSummary)=rownames(modulePreservIndVsConsensus$preservation$Z$ref.Consensus$inColumnsAlsoPresentIn.HSCC_H)
# colnames(preservSummary)=c("moduleSize", "Cons preserv in HSCC_H", "Cons preserv in HSCC_L")
# 
# write.csv(preservSummary, file="resultsCoexpr/coexprModulesPreserv.csv")
########################################################################################################
########################################################################################################
# save the results below for use with enrinchR
try(dir.create("resultsCoexpr_CEA/moduleGeneList"), silent = T)
coexprConnConsensus=intramodularConnectivity(adjCoexpr,  colorsCoexpr, scaleByMax=T)
#totalScaledConnectivity=coexprConnConsensus[,"kTotal"]/max(coexprConnConsensus[,"kTotal"])

coexprConnHigh=intramodularConnectivity(adj_HSCC_H,  colorsCoexpr, scaleByMax=T)
coexprConnLow=intramodularConnectivity(adj_HSCC_L,  colorsCoexpr, scaleByMax=T)

coexprResultsTable=cbind(colorsCoexpr, round(coexprConnConsensus[,"kWithin"],3), round(coexprConnHigh[,"kWithin"],3), round(coexprConnLow[,"kWithin"],3))
colnames(coexprResultsTable)=c("module", "consensus conn", "High Conn", "Low Conn")
rownames(coexprResultsTable)=geneNames

for (module in modulesCoexpr){
  print(module)
  currModuleInfo=cbind(rownames(coexprConnConsensus)[colorsCoexpr==module],round(coexprConnConsensus[colorsCoexpr==module,"kWithin"],2))
  colnames(currModuleInfo)=c("gene name", "module connectivity")
  write.csv(currModuleInfo, file=paste("resultsCoexpr_CEA/moduleGeneList/module_", module, ".csv", sep=""), row.names=F, col.names=F)  
}
#############################################################################
# GO annotations

load("/home/dan/workDir/networkAnalysis/data/transcriptInfoMouse.RData")
annotateMouseModulesGO(colorsCoexpr, transcriptInfoMouse, type="Coexpr_CEA")
##############################################################################
# record differential expression
results_highCPMgenes=read.csv("resultsCoexpr_CEA/resultsDEDV_highCPM.csv")
rownames(results_highCPMgenes)=results_highCPMgenes[,1]

intersect(geneNames, rownames(results_highCPMgenes))
setdiff(geneNames, rownames(results_highCPMgenes))

results_DEDVNet=results_highCPMgenes[geneNames,]


##############################################################################

##################################################################################33

# evaluate changes in edge strength
rawAdj1=adjacency(t(HSCC_H), power=1)
rawAdj2=adjacency(t(HSCC_L), power=1)

diffEdgesCEA = diffEdges(rawAdj1, rawAdj2, n1=dim(HSCC_H)[2], n2=dim(HSCC_L)[2], pThreshold=0.01, adjThreshold=0.5, nCores=7)
  
save(diffEdgesCEA, file="resultsCoexpr_CEA/diffEdgesCEA.RData")
load("resultsCoexpr_CEA/diffEdgesCEA.RData")

totalEdges=(length(geneNames))^2
affectedEdges=sum(diffEdgesCEA)
edgeChangeRate=affectedEdges/totalEdges

geneChangeEdgeCount=rowSums(diffEdgesCEA)
names(geneChangeEdgeCount)=geneNames

pValuesEdgeChange=rep(1,length(geneNames))
names(pValuesEdgeChange)=geneNames

for (gene in geneNames){
  pValuesEdgeChange[gene]=binom.test(x=geneChangeEdgeCount[gene], n=length(geneNames), p=edgeChangeRate, alternative  ="g")$p.value
}

#####################################################################################

pValues=pValuesEdgeChange
names(pValues)=geneNames

sortIndexes=sort.int(pValues, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pValues)
summary(adjustedResults)

sortedAdjustedPvals_DW=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals_DW)=sortedGeneNames

#####################################################################################

summaryResultsCEA=cbind(colorsCoexpr, results_DEDVNet[,c(1:3, 5:12)], round(coexprConnConsensus[,"kWithin"],3), round(coexprConnLow[,"kWithin"],3), round(coexprConnHigh[,"kWithin"],3),round(pValuesEdgeChange,3), round(sortedAdjustedPvals_DW[geneNames],3),geneChangeEdgeCount )
colnames(summaryResultsCEA)
colnames(summaryResultsCEA)=c("module","geneName", "logFC", "logCPM", "mean.counts.L","mean.counts.H","p.val.DE", "adj.p.DE" , "sd.counts.L","sd.counts.H","p.val.DV", "adj.p.DV","modular conn Consensus","modular conn Low", "modular conn High","p.val.DW", "adj.p.DW", "changed edge count")


#sanity checks
plot(summaryResultsCEA[,"logFC"], -log10(summaryResultsCEA[,"p.val.DE"]))
plot(summaryResultsCEA[,"changed edge count"], -log10(summaryResultsCEA[,"p.val.DW"]))
plot(summaryResultsCEA[,"changed edge count"], -log10(summaryResultsCEA[,"adj.p.DW"]))

write.csv(summaryResultsCEA, file="resultsCoexpr_CEA/tableResultsCEA.csv")


#####################################################################################
deGenes=geneNames[summaryResultsCEA[,"p.val.DE"] < 0.01]
dvGenes=geneNames[summaryResultsCEA[,"p.val.DV"] < 0.01]
dwGenes=geneNames[summaryResultsCEA[,"p.val.DW"] < 0.01]

modulesEnrichDE=moduleEnrichment(colorsCoexpr, deGenes)
affectedModulesDE=names(modulesEnrichDE)[modulesEnrichDE<(0.05/length(modulesCoexpr))]

modulesDE=""
for (i in 1:length(affectedModulesDE)){
  modulesDE=paste(modulesDE, affectedModulesDE[i], sep=" ")
}

fileConnSummary<-file("resultsCoexpr_CEA/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Modules affected by DE expression changes: ", modulesDE, sep=' '), fileConnSummary)
close(fileConnSummary)

############################################################################
modulesEnrichDV=moduleEnrichment(colorsCoexpr, dvGenes)
affectedModulesDV=names(modulesEnrichDV)[modulesEnrichDV<(0.05/length(modulesCoexpr))]

fileConnSummary<-file("resultsCoexpr_CEA/SummaryResultsCoexpr.txt",  open="at")

modulesDV=""
for (i in 1:length(affectedModulesDV)){
  modulesDV=paste(modulesDV, affectedModulesDV[i], sep=" ")
}

writeLines(paste("Modules affected by DV changes: ", modulesDV, sep=' '), fileConnSummary)
close(fileConnSummary)

############################################################################
modulesEnrichDW=moduleEnrichment(colorsCoexpr, dwGenes)
affectedModulesDW=names(modulesEnrichDW)[modulesEnrichDW<(0.05/length(modulesCoexpr))]

fileConnSummary<-file("resultsCoexpr_CEA/SummaryResultsCoexpr.txt",  open="at")
modulesDW=""
for (i in 1:length(affectedModulesDW)){
  modulesDW=paste(modulesDW, affectedModulesDW[i], sep=" ")
}

writeLines(paste("Modules affected by DW changes: ", modulesDW, sep=' '), fileConnSummary)
close(fileConnSummary)

############################################################################


# 
# names(colorsCoexpr)=geneNames
# edgesModuleEnrich=moduleEnrichment(colorsCoexpr, genesChangedEdges)
# affectedModules=names(edgesModuleEnrich[edgesModuleEnrich<(0.05/length(modulesCoexpr))])
# 
# fileConnSummary<-file("resultsCoexpr/SummaryResultsCoexpr.txt",  open="at")
# writeLines(paste("Modules affected by edge changes  ", affectedModules, sep=','), fileConnSummary)
# close(fileConnSummary)
# 
# 
# ##############################################################################
# ##############################################################################
# #calculate overlap with Mulligan 2006
# 
# MulliganTable=read.csv("data/SITable2_MulliganPonamarev_2006.csv", skip=1, header=T)
# MulliganGenes=as.character(MulliganTable[,1])
# 
# genesShellPval=summaryResultsShell[(summaryResultsShell[,"PValue"] <0.01 | summaryResultsShell[,"pvalVar"]<0.01 | summaryResultsShell[,"pValuesEdgeChange"]<0.01),"Gene.Names"]
# genesShellFdr=summaryResultsShell[(summaryResultsShell[,"fdrDE"] <0.01 | summaryResultsShell[,"fdrVar"]<0.01 | summaryResultsShell[,"fdrEdgesChange"]<0.01),"Gene.Names"]
# 
# intersect(MulliganGenes,genesShellPval)
# length(intersect(MulliganGenes,genesShellFdr))
# # length(intersect(MulliganGenes,genesShellFdr))
# # [1] 252
# fisher.test(summaryResultsShell[,"Gene.Names"] %in% genesShellPval, summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes, alternative = "g")
# fisher.test(summaryResultsShell[,"Gene.Names"] %in% genesShellFdr, summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes, alternative = "g")
# 
# MulliganDE=summaryResultsShell[,"Gene.Names"]  %in% MulliganGenes
# 
# summaryResultsShell=cbind(summaryResultsShell, MulliganDE)
# write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")
# 
# ##############################################################################
# hubsResultsShellCoexpr=read.csv("resultsCoexpr/hubsResultsShell.csv")[,"Gene.Names"]
# hubsResultsShellCoSplicEx=read.csv("resultsCoSplicEx/hubsResultsShell.csv")[,"Gene.Names"]
# 
# intersect(hubsResultsShellCoexpr, hubsResultsShellCoSplicEx)
# 
# ##############################################################################
# # summaryResultsShell=read.csv("resultsCoexpr/summaryResultsCoexprShell.csv")
# # 
# # summaryResultsShell=cbind(summaryResultsShell[,2:4], coexprConnHigh[,"kWithin"], coexprConnLow[,"kWithin"], summaryResultsShell[,6:18])
# # summaryResultsShell[,3:18]=round(summaryResultsShell[,3:18], 3)
# # colnames(summaryResultsShell)[1]="Gene Names"
# # colnames(summaryResultsShell)[4:5]=c("High Conn", "Low Conn")
# # write.csv(summaryResultsShell, file="resultsCoexpr/summaryResultsShell.csv")
# 
# ##########################################################################################################
# hubsResultsShell=read.csv("resultsCoexpr/hubsResultsShell.csv")
# 
# colorsHubsAffected=rep("white", length(geneNames))
# names(colorsHubsAffected)=geneNames
# colorsHubsAffected[hubsResultsShell[,"Gene.Names"]]="black"
# 
# # GO annotations
# 
# load("/home/dan/workDir/HDID2/data/transcriptInfoMouse.RData")
# annotateMouseModulesGO(colorsHubsAffected, transcriptInfoMouse, type="CoexprHubsAffected")
# 
