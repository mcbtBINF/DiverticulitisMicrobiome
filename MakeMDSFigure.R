#MakeMDSFigure
rm(list=ls())

#setwd("C:\\topeOneAtATime\\merged")
if (.Platform$OS.type == "windows") {
  baseDir <- c("/Google Drive/Datasets/Keku/Diverticulitis/")
} else {
  baseDir <- c("/Users/mbrown67/Google Drive/Datasets/Keku/Diverticulitis/")  
}

processedDir <- paste0(baseDir, "processed/")
analysisDir <- paste0(baseDir, "analysis/")

setwd(analysisDir)
#setwd(processedDir)
taxaLevels <- c("phylum", "class", "order", "family", "genus")

#taxa <- "genus"

for(taxa in taxaLevels ) 
{
  #taxa <- "genus"
myT <- read.table(paste0("mds_",taxa,"PlusMetadata.txt"), sep = "\t", header = TRUE)
myT <- myT[myT$numberSequencesPerSample >= 1000,]
#myT <- myT[myT$read == 1 &  myT$numberSequencesPerSample >= 1000 & !is.na(myT$caseControl)
#            & (myT$caseControl == 0 | myT$caseControl == 1 ), ]
myTeigenValues <- read.table(paste0("eigenValues_", taxa, "PlusMetadata.txt"), sep = "\t", header = TRUE)

myColors = ifelse(myT$caseControl == 1 , "red" , "black")

plot( myT$MDS1, myT$MDS2, col = myColors, pch = 16, cex = 1.2, 
      xlab=paste("MDS 1 (", format(myTeigenValues[1, 1], digits = 4), "%)", sep = ""), ylab = paste("MDS 2 (", format(myTeigenValues[2, 1], digits = 4), "%)", sep = ""))
print(wilcox.test(myT$MDS1[myT$caseControl == 0 ], myT$MDS1[myT$caseControl == 1 ]))$p.value
print(wilcox.test(myT$MDS2[myT$caseControl == 0 ], myT$MDS2[myT$caseControl == 1 ]))$p.value
}