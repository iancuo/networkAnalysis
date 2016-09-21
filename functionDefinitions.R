
###############################################################################################################################################

plotNetConstruction = function (sft)

{
sizeGrWindow(15, 9)
par(mfrow = c(2,3));
par(mar=c(5, 5, 3, 3))
par(lab=c(5, 7, 3))

cex1 = 1.5;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Fit",
     main = paste("Scale independence"), font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y.at <- c(seq(from = 0, to=1, by=.1))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)


abline(h=0.75,col="black")

mtext("A", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",type="l",
     main = paste("Mean Connectivity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,5]))

y_step=(y_range[2]-y_range[1])/7

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("B", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Density as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,8],
     xlab="Soft Threshold (power)",ylab="Density",type="l",
     main = paste("Density"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,8]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("C", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Centralization as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,9]))

plot(sft$fitIndices[,1], sft$fitIndices[,9],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Centralization"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");


x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,9]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("D", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)
############################################################################

# Heterogeneity as a function of the soft-thresholding power
y_range=c(0,max(sft$fitIndices[,10]))

plot(sft$fitIndices[,1], sft$fitIndices[,10],
     xlab="Soft Threshold (power)",ylab="Centralization",type="l", ylim=y_range,
     main = paste("Heterogeneity"),lwd=2, font.lab=2, cex.lab=1.2, xaxt = "n", yaxt = "n", col="red");

x.at <- powers
axis(1, at = x.at, labels = as.character(powers), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

y_range=c(0,max(sft$fitIndices[,10]))

y_step=(y_range[2]-y_range[1])/5

y.at <- c(seq(from = y_range[1], to=y_range[2], by=y_step))
axis(2, at = y.at, labels = format(y.at,digits=1), cex.lab=1.5, font=2, lwd=2, lwd.ticks=2)

mtext("E", NORTH<-3, at=-2, line=0.25, cex=2.5, font=2)

}

####################################################################################

annotateMouseModulesGO = function (colorsCoexpr, transcriptInfoMouse, type)
  # this assumes the gene identifiers are gene symbols, whichr are the names of the colorsCoexpr
  
{
  library("GOstats")
  
  dirName=paste("results", type, sep="")
  try(dir.create(dirName), silent = F)
  
  
  fileName=paste("results", type,"/modules", type, "GO.csv", sep="")
  modules=names(table(colorsCoexpr))
  geneNames=names(colorsCoexpr)    
  networkTranscriptInfo= transcriptInfoMouse[which(transcriptInfoMouse[,"external_gene_name"] %in% geneNames),]
  
  pValue=0.01/length(modules)
  univ= unique(as.character(networkTranscriptInfo[,"entrezgene"]))
  
  univ=unique(univ)
  
  for (module in modules){
    print(module)
    geneNamesModule=geneNames[colorsCoexpr==module]
    
    moduleEntrez = unique(as.character(networkTranscriptInfo[which(networkTranscriptInfo[,"external_gene_name"] %in% geneNamesModule),"entrezgene"]))
    BPparams <- new("GOHyperGParams", geneIds = moduleEntrez, universeGeneIds = univ,  annotation = "org.Mm.eg.db",ontology = "BP", pvalueCutoff = pValue, conditional = T, testDirection = "over")
        
    BPres <- hyperGTest(BPparams)
    res=summary(BPres)
    print(res)
    res=cbind(rep(paste(module,"module BP"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    
    CCparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "CC", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    CCres <- hyperGTest(CCparams)
    res=summary(CCres)
    res=cbind(rep.int(paste(module,"module CC"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    print(res)
    
    MFparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "MF", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    MFres <- hyperGTest(MFparams)
    res=summary(MFres)
    res=cbind(rep.int(paste(module,"module MF"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    print(res)
    
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
       
  }
    
}

####################################################################################

annotateMouseGeneGroupGO = function (affectedGenes,networkGenes, transcriptInfoMouse, type)
  # this assumes the gene identifiers are gene symbols, whichr are the names of the colorsCoexpr
  
{
  library("GOstats")
  
  dirName=paste("results", type, sep="")
  try(dir.create(dirName), silent = F)
  
  
  fileName=paste("results", type,"/affectedGenes", type, "GO.csv", sep="")
  networkTranscriptInfo= transcriptInfoMouse[which(transcriptInfoMouse[,"external_gene_name"] %in% networkGenes),]
  
  pValue=0.01
  univ= unique(as.character(networkTranscriptInfo[,"entrezgene"]))
  
  univ=unique(univ)
  module="affected genes"

    moduleEntrez = unique(as.character(networkTranscriptInfo[which(networkTranscriptInfo[,"external_gene_name"] %in% affectedGenes),"entrezgene"]))
    BPparams <- new("GOHyperGParams", geneIds = moduleEntrez, universeGeneIds = univ,  annotation = "org.Mm.eg.db",ontology = "BP", pvalueCutoff = pValue, conditional = T, testDirection = "over")
    
    BPres <- hyperGTest(BPparams)
    res=summary(BPres)
    print(res)
    res=cbind(rep(paste(module,"module BP"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    
    CCparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "CC", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    CCres <- hyperGTest(CCparams)
    res=summary(CCres)
    res=cbind(rep.int(paste(module,"module CC"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    
    print(res)
    
    MFparams <- new("GOHyperGParams", geneIds = moduleEntrez,
                    universeGeneIds = univ, annotation = "org.Mm.eg.db", 
                    ontology = "MF", pvalueCutoff = pValue, conditional = F,
                    testDirection = "over")
    
    MFres <- hyperGTest(MFparams)
    res=summary(MFres)
    res=cbind(rep.int(paste(module,"module MF"), dim(res)[1]),res)
    colnames(res)[1]="module and GO term"
    print(res)
    
    write.table(res, fileName, append=TRUE, sep =",", quote=F, row.names=F)
    

  
}



##########################################################################################
diffVarSplicing = function (canberraListSelected, phenotypeVector, nPerm, nCores)
{
   
  library(parallel)
  library(vegan)
  geneNames=names(canberraListSelected)
  samples=names(phenotypeVector)
  

  canberraListMatrices=foreach(gene=geneNames,.inorder=T,.verbose = F) %dopar% {
    
    as.matrix(canberraListSelected[[gene]])[samples, samples]

  } 
  names(canberraListMatrices)=geneNames
  
  mrppResults=vector(mode = "list", length(geneNames))
  names(mrppResults)=geneNames
  
  mrppA=rep(0, length(geneNames))
  names(mrppA)=geneNames
  mrppPvalues=mrppA
  
  
  for (gene in geneNames) {
    print(which(geneNames==gene))
    if (sum(is.na(canberraListMatrices[[gene]]))==0){
      mrppResults[[gene]]=mrpp(canberraListMatrices[[gene]], grouping=phenotypeVector, permutations=nPerm, parallel=nCores)
      mrppA[gene]=mrppResults[[gene]]$A
      mrppPvalues[gene]=mrppResults[[gene]]$Pvalue
      
    }
  }
    
  rm(mrppResults)
  sortIndexes=sort.int(mrppPvalues, decreasing = F, index.return=T)$ix
  sortedGeneNames=geneNames[sortIndexes]
  
  adjustedResults<-SGoF(u=mrppPvalues)
  
  
  summary(adjustedResults)
  
  sortedPvals=mrppPvalues[sortedGeneNames]
  sortedA=mrppA[sortedGeneNames]
  sortedAdjustedPvals=adjustedResults$Adjusted.pvalues
  names(sortedAdjustedPvals)=names(sortedPvals)
  
  
  returnValue=cbind(sortedA, sortedPvals, sortedAdjustedPvals)
  return(returnValue)
}

##########################################################################################


moduleEnrichment = function (colors, signifGenes)
  
{ 
  geneNames=names(colors)
  modules=names(table(colors))
  modulePvaluesEnrich=rep(1,length(modules))
  names(modulePvaluesEnrich)=modules
  for (module in modules){
    moduleGenes=geneNames[colors==module]
    modulePvaluesEnrich[module]=fisher.test(geneNames %in% moduleGenes, geneNames %in% signifGenes, alternative = "g")$p.value
  }
  return(modulePvaluesEnrich)
}


#############################################################################################
##########################################################################################


diffEdges = function (rawAdj1, rawAdj2, n1, n2, pThreshold=0.01, adjThreshold=0.5, nCores){
  library(psych)
  
  options(cores=nCores)
  geneNames=rownames(rawAdj1)
  nGenes=length(geneNames)
  
  diffEdgesMatrix <- foreach (gene=geneNames, .inorder=T, .verbose = T, .combine=cbind) %dopar% {
    
 
    vectorEdges=rep(1, nGenes)
    names(vectorEdges)=geneNames
    
    vector1=rawAdj1[gene,]
    vector2=rawAdj2[gene,]
    idxsVector=abs(vector1 -vector2) > adjThreshold
    
    for (gene2 in geneNames[idxsVector==T]){
      vectorEdges[gene2]=r.test(n=n1, r12=rawAdj1[gene,gene2], r34=rawAdj2[gene,gene2], n2=n2)$p
    }
    
    vectorEdges<pThreshold 
  }
  
  return(diffEdgesMatrix)
}

