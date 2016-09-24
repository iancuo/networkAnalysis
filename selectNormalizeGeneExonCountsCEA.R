      
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(sgof)


getDoParWorkers()
options(cores=14)
getDoParWorkers()

# linux does not work with box so this code has a local data directory

setwd("/home/dan/workDir/networkAnalysis")

source("/home/dan/workDir/networkAnalysis/functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = F)
try( dir.create("figuresCoexpr"), silent = F)
try(dir.create("resultsCoSplicEx"), silent = F)
try( dir.create("figuresCoSplicEx"), silent = F)

# read raw data - can be improved by using read.table from original .txt file
geneReadsRaw=read.table("RNASeq019 HSCC Alex/CEA/RNASeq019_CEA_mm10_gene_reads_not_normalized.txt")
geneNames=rownames(geneReadsRaw)

# read sample info - sample names need to be inspected and categories extracted differently for each dataset !!!!!!!!!
sampleKey=read.csv("RNASeq019 HSCC Alex/CEA/CEAsampleKey.csv", header=T)
sampleKey[,"CeA"]=paste("S", sampleKey[,"CeA"], sep = "")

# other projects might mean strsplit
# splitIDs=mapply(strsplit, sample_info, MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
# sampleInfo=unlist(lapply(splitIDs, "[[", 3))

samplesHigh=sampleKey[grep("H",sampleKey[,"Line"] ), "CeA"]
samplesLow=sampleKey[grep("L",sampleKey[,"Line"] ), "CeA"]

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
# results from sgof come out sorted so the original pvalues and geneNames need to be sorted

pValues=de.tgw$table$PValue
names(pValues)=geneNames

sortIndexes=sort.int(pValues, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pValues)
summary(adjustedResults)

sortedAdjustedPvals=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals)=sortedGeneNames

# some sanity checks
#plot(pValues[sortedGeneNames], adjustedResults$data) # should be straight line
#plot(pValues[sortedGeneNames], adjustedResults$data) # should be straight line

meanCounts=rowMeans(geneReadsRaw)
sdCounts=apply(geneReadsRaw,1, sd) 
cvCounts=sdCounts/meanCounts

resultsDEtotal=cbind(de.tgw$table[sortedGeneNames,], meanCounts[sortedGeneNames], cvCounts[sortedGeneNames] )
colnames(resultsDEtotal)=c(colnames(de.tgw$table), c("mean counts", "CV counts"))

#this will be collected in Supplemental Table 1
write.csv(resultsDEtotal, file="resultsCoexpr/resultsDEtotal.csv")


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


# select genes with logCPM > 0 (equivalent to CPM>1) and with high CV (top 75%) for further analysis

geneNamesHighCPM=geneNames[resultsDEtotal[,"logCPM"]>0]
geneCountsHighCPM=normalizedGeneCountsUQ[geneNamesHighCPM,]
cvCountsHighCPM=cvCounts[geneNamesHighCPM]

quantileCV=quantile(cvCountsHighCPM, na.rm=T)  
geneNamesHighCV=geneNamesHighCPM[cvCountsHighCPM>quantileCV[2]]
geneReadsHighCV=geneCountsHighCPM[cvCountsHighCPM>quantileCV[2],]


#############################################################################################################
connCoexpr=rowSums(geneReadsHighCV)
hist(connCoexpr,100)
quantile(connCoexpr, na.rm=T)


plot(cvCountsHighCPM[geneNamesHighCV], connCoexpr, ylim=c(0,500000), xlab="CV gene counts", ylab="Coexpr connectivity", cex.lab=1.25, font.lab=2)

cor(cvCountsHighCPM[geneNamesHighCV], connCoexpr, use = "pairwise.complete.obs")

quantileConnGenes=quantile(connCoSplicEx, probs = seq(0, 1, 0.25))  
quantileCVdist=quantile(cvDist, probs = seq(0, 1, 0.25), na.rm=T)  

plot(cvDist, connCoSplicEx, xlim=c(0,1), xlab="CV Canberra distances", ylab="CoSplicEx connectivity", cex.lab=1.25, font.lab=2)
hist(cvDist, 100)
hist(connCoSplicEx, 50)
sum(connCoSplicEx > 2)
sum(connCoSplicEx > 1.5)


geneNamesHighCoSplicExConn=names(canberraListExons)[connCoSplicEx>quantileConnExons[6]]
################################################################################################3



###################################################################################
exonCounts=read.table("RNASeq019 HSCC Alex/CEA/RNASeq019_CEA_mm10_exon_reads_not_normalized.txt")


splitIDs=mapply(strsplit, rownames(exonCounts), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))
exon_start=unlist(lapply(splitIDs, "[[", 2))

#sanity check

colnames(exonCounts)==colnames(geneReadsHighCV)

sum(exonCounts[exonGeneName=="Drd2",])
sum(geneReadsHighCV["Drd2",])

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
save(canberraListExons, file="resultsCoexpr/canberraListExonsCEA.RData")
# load("data/canberraListExons.RData")
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
# compute CV of pairwise distances

meanDist=rowMeans(t(distData))
sdDist=apply(distData,2, sd) 
cvDist=sdDist/meanDist
#############################################################################################################
adjCoSplicEx_large=adjacency(distData,power=6)
save(adjCoSplicEx_large, file="resultsCoexpr/adjCoSplicEx_large.RData")

#load("data/adjCoSplicEx_large.RData")
#just in case ...
adjCoSplicEx_large[is.na(adjCoSplicEx_large)]=0
diag(adjCoSplicEx_large)=1
colnames(adjCoSplicEx_large)=rownames(adjCoSplicEx_large)
connCoSplicEx=rowSums(adjCoSplicEx_large)


#############################################################################################################


quantileConnExons=quantile(connCoSplicEx, probs = seq(0, 1, 0.25))  
quantileCVdist=quantile(cvDist, probs = seq(0, 1, 0.25), na.rm=T)  

plot(cvDist, connCoSplicEx, xlim=c(0,1), xlab="CV Canberra distances", ylab="CoSplicEx connectivity", cex.lab=1.25, font.lab=2)
hist(cvDist, 100)
hist(connCoSplicEx, 50)
sum(connCoSplicEx > 2)
sum(connCoSplicEx > 1.5)


geneNamesHighCoSplicExConn=names(canberraListExons)[connCoSplicEx>quantileConnExons[6]]
################################################################################################3
selectedGeneCounts=geneReadsHighConn

canberraListSelected=canberraListExons[geneNamesHighCoSplicExConn]
adjCoSplicEx=adjCoSplicEx_large[geneNamesHighCoSplicExConn,geneNamesHighCoSplicExConn]

exonGeneNameSelected=geneNamesHighCoSplicExConn
selectedExonCounts=normExonCounts[which(exonGeneName %in% geneNamesHighCoSplicExConn),]

########################################################################################################

save(selectedGeneCounts, canberraListSelected,adjCoSplicEx,selectedExonCounts, exonGeneNameSelected, groupSelection, samplesHigh, samplesLow, file="data/selectedData.RData")
load("data/selectedData.RData")



############################################################################################
