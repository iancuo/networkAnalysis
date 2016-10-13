library(matrixStats)      
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(sgof)
library(multtest)

getDoParWorkers()
options(cores=14)
getDoParWorkers()

# linux does not work with box so this code has a local data directory

setwd("/home/dan/workDir/networkAnalysis")

source("/home/dan/workDir/networkAnalysis/functionDefinitions.R")
try(dir.create("resultsCoexpr_SH"), silent = F)
try( dir.create("figuresCoexpr_SH"), silent = F)
try(dir.create("resultsCoSplicEx_SH"), silent = F)
try( dir.create("figuresCoSplicEx_SH"), silent = F)

# read raw data - can be improved by using read.table from original .txt file
geneReadsRaw=as.matrix(read.table("RNASeq019 HSCC Alex/SH/RNASeq019_SH_mm10_gene_reads_not_normalized.txt"))

geneNames=rownames(geneReadsRaw)

# read sample info - sample names need to be inspected and categories extracted differently for each dataset !!!!!!!!!
sampleKey=read.csv("RNASeq019 HSCC Alex/SH/SHsampleKey.csv", header=T)
sampleKey[,"Shell"]=paste("S", sampleKey[,"Shell"], sep = "")

# other projects might mean strsplit
# splitIDs=mapply(strsplit, sample_info, MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
# sampleInfo=unlist(lapply(splitIDs, "[[", 3))

samplesHigh=sampleKey[grep("H",sampleKey[,"Line"] ), "Shell"]
samplesLow=sampleKey[grep("L",sampleKey[,"Line"] ), "Shell"]

##################################################################################################################
# divide the data in different groups and perform DE with edgeR
HSCC_H=geneReadsRaw[,samplesHigh ]
HSCC_L=geneReadsRaw[,samplesLow]

groupSelection=c(rep("HSCC_H",dim(HSCC_H)[2]),rep("HSCC_L",dim(HSCC_L)[2]))
groupSelection =factor(groupSelection)

d=DGEList(counts= cbind(HSCC_H, HSCC_L), group= groupSelection)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 
########################################################################################

# use sgof package for multiple comparison correction
# results from sgof come out sorted but un-named (!!!!) so the original pvalues and geneNames need to be sorted

pValues_DE=de.tgw$table$PValue
names(pValues_DE)=geneNames

sortIndexes=sort.int(pValues_DE, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pValues_DE)
summary(adjustedResults)

sortedAdjustedPvals_DE=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals_DE)=sortedGeneNames

fileConnSummary<-file("resultsCoexpr_SH/SummaryResultsCoexpr.txt",  open="wt")

writeLines(paste("Number of genes with >1 CPM that are DE at FDR=0.05: ", sum(sortedAdjustedPvals_DE<0.05), sep=' '), fileConnSummary)
close(fileConnSummary)

# some sanity checks
#plot(pValues[sortedGeneNames], adjustedResults$data) # should be straight line
#plot(pValues[sortedGeneNames], adjustedResults$data) # should be straight line

##############################################################################
###################################################################################

# calculate edgeR normalization factors and normalize the data
UQnormFactors=calcNormFactors(geneReadsRaw, method=c("upperquartile"))

effectiveLibrarySizes= UQnormFactors*colSums(geneReadsRaw)
meanEffLibSize=mean(effectiveLibrarySizes)
countNormFactor= meanEffLibSize/effectiveLibrarySizes

normalizedGeneCountsUQ=0* geneReadsRaw

for (sample in 1:dim(normalizedGeneCountsUQ)[2]){	
	normalizedGeneCountsUQ[,sample]= geneReadsRaw[, sample]* countNormFactor[sample]	
}

# verify that smaller libraries are multiplied by bigger normalization factors
plot(colSums(geneReadsRaw), countNormFactor, xlab="Un-normalized library sizes", ylab="Normalization factors")

# plot original library sizes versus normalized
#range of data much smaller in normalized data
xylim=c(min(colSums(geneReadsRaw), colSums(normalizedGeneCountsUQ)), max(colSums(geneReadsRaw), colSums(normalizedGeneCountsUQ)))
plot(colSums(geneReadsRaw), colSums(normalizedGeneCountsUQ), xlim=xylim, ylim=xylim)

# interquartile range of normalized libraries should me much smaller
IQR(colSums(geneReadsRaw))
IQR(colSums(normalizedGeneCountsUQ))

###################################################################################

###################################################################################
# select genes with logCPM > 0 (equivalent to CPM>1) 

geneNamesHighCPM=geneNames[de.tgw$table[,"logCPM"]>0]
geneCountsHighCPM=normalizedGeneCountsUQ[geneNamesHighCPM,]

meanCounts_H=rowMeans(HSCC_H[geneNamesHighCPM,])
meanCounts_L=rowMeans(HSCC_L[geneNamesHighCPM,])

sdCounts_H=apply(HSCC_H[geneNamesHighCPM,],1, sd) 
sdCounts_L=apply(HSCC_L[geneNamesHighCPM,],1, sd) 


##############################################################################
# find differentially variable genes

pvalVar=rep(1, length(geneNamesHighCPM))
names(pvalVar)=geneNamesHighCPM

for (gene in geneNamesHighCPM){
  pvalVar[gene]=var.test(x=HSCC_H[gene,], y=HSCC_L[gene,])$p.value
}

pvalVar[is.na(pvalVar)]=1

pValues=pvalVar
names(pValues)=geneNamesHighCPM

sortIndexes=sort.int(pValues, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNamesHighCPM[sortIndexes]

adjustedResults<-SGoF(u=pValues)
summary(adjustedResults)

sortedAdjustedPvals_DV=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals_DV)=sortedGeneNames

fileConnSummary<-file("resultsCoexpr_SH/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Number of genes with >1 CPM that are DV at FDR=0.05: ", sum(sortedAdjustedPvals_DV<0.05), sep=' '), fileConnSummary)
close(fileConnSummary)

geneNamesDE=sortedGeneNames[sortedAdjustedPvals_DE < 0.05]
geneNamesDV=sortedGeneNames[sortedAdjustedPvals_DV < 0.05]

write.csv(geneNamesDE, file="resultsCoexpr_SH/geneNamesDE.csv")
write.csv(geneNamesDV, file="resultsCoexpr_SH/geneNamesDV.csv")

setdiff(geneNamesHighCPM, geneNames)

###########################################################################33
results_highCPMgenes=cbind(de.tgw$table[geneNamesHighCPM,], meanCounts_L[geneNamesHighCPM], meanCounts_H[geneNamesHighCPM], pValues_DE[geneNamesHighCPM], sortedAdjustedPvals_DE[geneNamesHighCPM],sdCounts_L[geneNamesHighCPM],sdCounts_H[geneNamesHighCPM],  pvalVar, sortedAdjustedPvals_DV[geneNamesHighCPM] )
results_highCPMgenes=round(results_highCPMgenes,3)
colnames(results_highCPMgenes)=c(colnames(de.tgw$table), c("mean counts L", "mean counts H", " p val DE", " adj p DE", "sd L", "sd H", "p val DV", "adj p val DV"))
rownames(results_highCPMgenes)=geneNamesHighCPM
#this will be collected in Supplemental Table 1
write.csv(results_highCPMgenes, file="resultsCoexpr_SH/resultsDEDV_highCPM.csv")

#######################################################################################3
# possibly swithch to bicor correlation in the future
#adjCoexpr=adjacency(t(geneCountsHighCPM), corFnc = "bicor", type="unsigned", power=6)

adjCoexpr=adjacency(t(geneCountsHighCPM),  type="unsigned", power=6)
adjCoexpr[is.na(adjCoexpr)]=0

#sanity check, should be 0
sum(adjCoexpr<0, na.rm=T)
connCoexpr=rowSums(adjCoexpr, na.rm = T)-1
connCoexpr_WGCNA=softConnectivity(t(geneCountsHighCPM), type="unsigned", power=6)
plot(connCoexpr, connCoexpr_WGCNA)

names(connCoexpr)=geneNamesHighCPM

connCoexpr=connCoexpr[connCoexpr > 0]

connCoexpr[is.na(connCoexpr)]=0
sortedConn=sort(connCoexpr, decreasing = T)
sum(connCoexpr < 0, na.rm=T)

totalConn=sum(sortedConn)
cumulativeConnFraction=0*sortedConn
for (i in 1:length(sortedConn)){
  cumulativeConnFraction[i]=sum(sortedConn[1:i])/totalConn
}


lastGene=min(which(cumulativeConnFraction > .9))

highConnGenes=names(sortedConn)[1:lastGene]

plot(cumulativeConnFraction, xlab="", ylab="")
title(xlab="Number of genes \n ranked by connectivity", cex.lab=1.25, font.lab=2)
title(ylab="Fraction total connectivity", cex.lab=1.25, font.lab=2)
abline(v=lastGene)
abline(h=0.9)
text(x=8100, y=0.8, labels=paste(lastGene, " genes", sep=""))
text(x=3100, y=0.95, labels="90% connectivity captured")
title(main="Selecting the network size\n starting from genes with > 1 CPM", cex.lab=1.5, font.lab=2)
# export figure
#############################################################################################################
adjCoexprHighConn=adjCoexpr[highConnGenes,highConnGenes]

selectedGeneCounts=normalizedGeneCountsUQ[highConnGenes,]
###################################################################################
exonCounts=read.table("RNASeq019 HSCC Alex/SH/RNASeq019_SH_mm10_exon_reads_not_normalized.txt")


splitIDs=mapply(strsplit, rownames(exonCounts), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))
exon_start=unlist(lapply(splitIDs, "[[", 2))

#sanity check

colnames(exonCounts)==colnames(geneReadsRaw)

# sanity check
sum(exonCounts[exonGeneName=="Drd2",])
sum(selectedGeneCounts["Drd2",])

normExonCounts=0* exonCounts
for (sample in 1:dim(exonCounts)[2]){
  normExonCounts[,sample]= exonCounts[, sample]* countNormFactor[sample]	
}
normExonCounts =round(normExonCounts)
rownames(normExonCounts)=rownames(exonCounts)

# select exons from genes with at least 1 CPM
exonCountsHighCounts=normExonCounts[which(exonGeneName %in% geneNamesHighCPM),]
exonGeneNamesHighCounts=exonGeneName[exonGeneName %in% geneNamesHighCPM]

# compute Canberra distances 
canberraListExons=foreach (geneName = geneNamesHighCPM, .inorder=T, .verbose = T) %dopar% {
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  # colnames(exonDistMatrix)=exonColnames
  # rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExons)=geneNamesHighCPM
save(canberraListExons, file="resultsCoSplicEx_SH/canberraListExonsSH.RData")
 load("resultsCoSplicEx_SH/canberraListExonsSH.RData")
########################################################################################################

nGenes=length(canberraListExons)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExons[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExons))
colnames(distData)=names(canberraListExons)

for(gene in names(canberraListExons)) {
  distData[,gene]=as.vector(as.dist(canberraListExons[[gene]]))
} 

##########################################################################################################
# compute connectivity in large CoSplicEx network

#######################################################################################3
adjCoSplicEx_raw=adjacency(distData,  type="unsigned", power=6)

adjCoSplicEx_raw[is.na(adjCoSplicEx_raw)]=0

#sanity check, should be 0
sum(adjCoSplicEx_raw<0, na.rm=T)
connCoSplicEx=rowSums(adjCoSplicEx_raw, na.rm = T)-1

names(connCoSplicEx)=geneNamesHighCPM

connCoSplicEx=connCoSplicEx[connCoexpr > 0]

connCoSplicEx[is.na(connCoSplicEx)]=0
sortedConn=sort(connCoSplicEx, decreasing = T)
sum(connCoSplicEx < 0, na.rm=T)

totalConn=sum(sortedConn)
cumulativeConnFraction=0*sortedConn
for (i in 1:length(sortedConn)){
  cumulativeConnFraction[i]=sum(sortedConn[1:i])/totalConn
}


lastGene=min(which(cumulativeConnFraction > .9))

highConnGenes=names(sortedConn)[1:lastGene]

plot(cumulativeConnFraction, xlab="", ylab="")
title(xlab="Number of genes \n ranked by connectivity", cex.lab=1.25, font.lab=2)
title(ylab="Fraction total connectivity", cex.lab=1.25, font.lab=2)
abline(v=lastGene)
abline(h=0.9)
text(x=6100, y=0.8, labels=paste(lastGene, " genes", sep=""))
text(x=3100, y=0.95, labels="90% connectivity captured")
title(main="Selecting the network size\n starting from genes with > 1 CPM", cex.lab=1.5, font.lab=2)
# export figure
#############################################################################################################
adjCoSplicEx=adjCoSplicEx_raw[highConnGenes,highConnGenes]


canberraListSelected=canberraListExons[highConnGenes]

exonGeneNameSelected=highConnGenes
selectedExonCounts=normExonCounts[which(exonGeneName %in% exonGeneNameSelected),]

########################################################################################################

save(selectedGeneCounts, canberraListSelected,adjCoSplicEx,selectedExonCounts, exonGeneNameSelected, groupSelection, samplesHigh, samplesLow, file="selectedData_SH.RData")
#load("data/selectedData.RData")



############################################################################################
