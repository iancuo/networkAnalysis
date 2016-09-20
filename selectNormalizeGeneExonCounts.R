      
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

setwd("/home/dan/workDir/networkAnalysis")

source("/home/dan/workDir/functionDefinitions.R")
try(dir.create("resultsCoexpr"), silent = F)
try( dir.create("figuresCoexpr"), silent = F)
try(dir.create("resultsCoSplicEx"), silent = F)
try( dir.create("figuresCoSplicEx"), silent = F)

# read raw data - can be improved by using read.table from original .txt file
geneReadsRaw=as.matrix(read.csv("data/RNA160225TP_gene_reads_not_normalized.csv", header=T, row.names=1))
geneNames=rownames(geneReadsRaw)
# read sample info - sample names need to be inspected and categories extracted differently for each dataset !!!!!!!!!
sample_info=colnames(geneReadsRaw)
splitIDs=mapply(strsplit, sample_info, MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
sampleInfo=unlist(lapply(splitIDs, "[[", 3))

samplesHigh=grep("H",sampleInfo )
samplesLow=grep("L",sampleInfo )

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
names(sortedAdjustedPvals)=names(sortedPvals)

resultsDEtotal=cbind(de.tgw$table, de.calls)

# select genes with logCPM > 0 (equivalent to CPM>1) for further analysis
resultsDEtotal=resultsDEtotal[resultsDEtotal[,"logCPM"]>0,]
write.csv(resultsDEtotal, file="resultsCoexpr/resultsDEtotal.csv")
# resultsDEtotal=read.csv("resultsCoexpr/resultsDEtotal.csv")
# rownames(resultsDEtotal)=resultsDEtotal[,1]

###################################################################################

# calculate edgeR normalization factors and normalize the data - use all data not just selected
UQnormFactors=calcNormFactors(geneReadsRaw, method=c("upperquartile"))

effectiveLibrarySizes= UQnormFactors*colSums(geneReadsRaw)
meanEffLibSize=mean(effectiveLibrarySizes)
countNormFactor= meanEffLibSize/effectiveLibrarySizes

normalizedGeneCountsUQ=0* geneReadsRaw

for (sample in 1:dim(normalizedGeneCountsUQ)[2]){	
	normalizedGeneCountsUQ[,sample]= geneReadsRaw[, sample]* countNormFactor[sample]	
}

geneReads=normalizedGeneCountsUQ[rownames(resultsDEtotal), ]
##################################################
# # find genes with connectivity in top 50%

CPMgenes=rownames(geneReads)[resultsDEtotal[,"logCPM"] > 0]

connCounts=softConnectivity(t(geneReads), power=6)

quantileConn=quantile(connCounts, seq(0, 1, 0.1))  
geneReadsHighConn=geneReads[connCounts>quantileConn[6],]
geneNamesHighConn=rownames(geneReadsHighConn)
###################################################################################
exonCounts=read.csv("data/RNA160225TP_exon_reads_not_normalized.csv", header=T, row.names=1)

splitIDs=mapply(strsplit, rownames(exonCounts), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))
exon_start=unlist(lapply(splitIDs, "[[", 2))

sum(exonCounts[exonGeneName=="Drd2",])
sum(geneReadsHighConn[geneNamesHighConn=="Drd2",])

exon_unique_id=mapply(paste, exonGeneName, exon_start, MoreArgs=list(sep="_"))
names(exonGeneName)=exon_unique_id

rownames(exonCounts)= exon_unique_id
sampleNames=as.vector(colnames(exonCounts))

normExonCounts=0* exonCounts
for (sample in 1:dim(exonCounts)[2]){
  normExonCounts[,sample]= exonCounts[, sample]* countNormFactor[sample]	
}
normExonCounts =round(normExonCounts)
rownames(normExonCounts)=rownames(exonCounts)


# select exons from genes with at least 1 CPM
exonCountsHighCounts=normExonCounts[which(exonGeneName %in% rownames(resultsDEtotal)),]
exonGeneNamesHighCounts=exonGeneName[exonGeneName %in% rownames(resultsDEtotal)]
  
canberraListExons=foreach (geneName = rownames(resultsDEtotal), .inorder=T, .verbose = T) %dopar% {
  #geneName=geneNames[i]
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
names(canberraListExons)=rownames(resultsDEtotal)

########################################################################################################


########################################################################################################

save(canberraListExons, file="data/canberraListExons.RData")
load("data/canberraListExons.RData")

nGenes=length(canberraListExons)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExons[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExons))
colnames(distData)=names(canberraListExons)

for(gene in names(canberraListExons)) {
  distData[,gene]=as.vector(as.dist(canberraListExons[[gene]]))
} 

adjCoSplicEx_large=adjacency(distData,power=6)


save(adjCoSplicEx_large, file="data/adjCoSplicEx_large.RData")

#load("data/adjCoSplicEx_large.RData")
#just in case ...
adjCoSplicEx_large[is.na(adjCoSplicEx_large)]=0
diag(adjCoSplicEx_large)=1
colnames(adjCoSplicEx_large)=rownames(adjCoSplicEx_large)


connCoSplicEx=rowSums(adjCoSplicEx_large)
quantileConnExons=quantile(connCoSplicEx, probs = seq(0, 1, 0.1))  

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
