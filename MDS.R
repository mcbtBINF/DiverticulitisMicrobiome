# Code for the MDS


rm(list = ls())
library("vegan")

#setwd("C:\\topeOneAtATime\\merged")
if (.Platform$OS.type == "windows") {
  baseDir <- c("/Google Drive/Datasets/Keku/Diverticulitis/")
} else {
  baseDir <- c("/Users/mbrown67/Google Drive/Datasets/Keku/Diverticulitis/")  
}

processedDir <- paste0(baseDir, "processed/")
analysisDir <- paste0(baseDir, "analysis/")

taxaLevels <- c("phylum", "class", "order", "family", "genus")

for(taxa in taxaLevels ) 
{
  setwd(processedDir)
  inFileName <- paste("pivoted_", taxa, "asColumnsLogNormalPlusMetadata.txt", sep = "")
  myT <-read.table(inFileName, header=TRUE, sep="\t")
  #numCols <- ncol(myT)
  #myColClasses <- c("character", rep("numeric", numCols - 1))
  #myT <-read.table(inFileName, header=TRUE, sep="\t",row.names=1, colClasses = myColClasses)
  myT <- myT[myT$read == 1,]
  myT <- myT[myT$caseControl != -1,]
  myPCOA <- capscale(myT[,(20:ncol(myT))]~1,distance="bray")
  
  setwd(analysisDir)
  MDSwithMetadata <- cbind(myT[,1:20], myPCOA$CA$u)
  write.table(MDSwithMetadata, sep = "\t", file = paste("mds_", taxa, "PlusMetadata.txt",sep = ""))
  write.table(myPCOA$CA$eig, file = paste("eigenValues_", taxa, "PlusMetadata.txt", sep = ""), sep = "\t")
}