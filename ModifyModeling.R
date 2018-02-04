# Adapted from A. Fodor's code 
# Available at: https://github.com/afodor/metagenomicsTools/blob/master/src/scripts/topeOneAtATime/metadataTests.txt

rm(list=ls())
library("Kendall")
library("vegan")
library("coin")

if (.Platform$OS.type == "windows") {
  baseDir <- c("/Google Drive/Datasets/Keku/Diverticulitis/")
} else {
  baseDir <- c("/Users/mbrown67/Google Drive/Datasets/Keku/Diverticulitis/")  
}

processedDir <- paste0(baseDir, "processed/")
analysisDir <- paste0(baseDir, "analysis/")


taxaLevels <- c("phylum", "class", "order", "family", "genus", "otu")
#taxaLevels <- c("phylum","class","order","family","genus","otu_qiime_cr","otu")
#taxaLevels <- c("genus")

for(taxa in taxaLevels ) {
  setwd(processedDir)
  
  inFileName <- paste("pivoted_", taxa, "asColumnsPlusMetadata.txt", sep = "")
  myT <-read.table(inFileName, header = TRUE, sep = "\t")
  myT <- myT[ myT$read == 1 &  myT$numberSequencesPerSample >= 1000 & !is.na(myT$caseControl), ]
  
  ## Dropping "blank" samples which correspond to caseControl value of -1
  myT <- myT[myT$caseControl != -1,]
  ifelse(taxa == "otu", startCol <- 19, startCol <- 20 )
  
  speciesRichness <- rarefy(myT[ , startCol:ncol(myT)], sample = min(myT$numberSequencesPerSample))

  inFileName <- paste("pivoted_", taxa, "asColumnsLogNormalPlusMetadata.txt", sep = "")
  myT <-read.table(inFileName,header = TRUE, sep = "\t")
  myT <- myT[ myT$read == 1 &  myT$numberSequencesPerSample >= 1000 & !is.na(myT$caseControl), ]
  
  ## Dropping "blank" samples which correspond to caseControl value of -1
  myT <- myT[myT$caseControl != -1, ]
  
  setwd(analysisDir)
  
  
  #myT$caseControl <- as.factor(myT$caseControl)
  
  sampleSizeCase <- vector()
  sampleSizeControl <- vector()
  pValuesWaist <- vector()
  rSquaredWaist <- vector()
  pValuesWaistKendall <- vector()
  effectSizeWaistKendall <- vector()
  pValuesWHR <- vector()
  rSquaredWHR <- vector()
  pValuesWHRKendall <- vector()
  effectSizeWHRKendall <- vector()
  pValuesTicsCount <- vector()
  rSquaredTicsCount <- vector()
  pValuesTicsCountKendall <- vector()
  effectSizeTicsCountKendall <- vector()
  pValuesWBO <- vector()
  rSquaredWBO <- vector()
  pValuesWBOKruskalWallis <- vector()
  effectSizeWBOKruskalWallis <- vector()
  ## Nonparametric vectors here
  pValuesCaseControl <- vector()
  rSquaredCaseControl <- vector()
  pValuesCaseControlWilcox <- vector()
  rSquaredCaseControlWilcox <- vector()
  pValuesAge <- vector()
  rSquaredAge <- vector()
  pValuesAgeKendall <- vector()
  effectSizeAgeKendall <- vector()
  pValuesSex <- vector()
  rSquaredSex <- vector()
  pValuesSexWilcox <- vector()
  rSquaredSexWilcox <- vector()
  #ticLocationPValue <- vector()
  pValuesticLocation <- vector()
  rSquaredticLocation <- vector()
  #ticLocationRSquared <- vector()
  pValuesTicLocationWilcox <- vector()
  rSquaredTicLocationWilcox <- vector()
  pValuesHemSize <- vector()
  rSquaredHemSize <- vector()
  pValuesHemSizeWilcox <- vector()
  rSquaredHemSizeWilcox <- vector()
  pValuesBMI <- vector()
  rSquaredBMI <- vector()
  pValuesBMIKendall <- vector()
  effectSizeBMIKendall <- vector()
  pValuesBMICAT <- vector()
  rSquaredBMICAT <- vector()
  pValuesBMICATKruskalWallis <- vector()
  effectSizeBMICATKruskalWallis <- vector()
  ## Nonparametric for BMICAT
  pValuesShannon <- vector()
  rSquaredShannon <- vector()
  rSquaredShannonWilcox <- vector()
  pValuesRichness <- vector()
  rSquaredRichness <- vector()
  rSquaredRichnessWilcox <- vector()
  rSquaredRichnessWilcox <- vector()
  meanCase <- vector()
  meanControl <- vector()
  names <- vector()
  #namesTicLocation <- vector()
  #namesHemSize <- vector()
  #print("Fails at Shannon Diversity")
  ShannonDiversity <- apply(myT[ , startCol:ncol(myT)], 1, diversity)
  #myT <- cbind(myT, ShannonDiversity, speciesRichness)

  pdf(paste("boxplots_", taxa, ".pdf", sep = ""))
  par(mfrow = c(2, 2))
  
  index <- 1

  for( i in c(startCol:(ncol(myT)))) 	
    if( sum( myT[ , i] > 0 ) > nrow(myT) / 4 ) 
    {	
      ## For easier use of wilcox_test
      bug <- myT[ , i]
      caseControl <- myT$caseControl
      sex <- myT$sex
      wbo <- myT$wbo
      bmi_CAT <- myT$bmi_CAT
      ticLocation <- myT$ticLocation
      #hemsize_s_ml <- myT$hemsize_s_ml
      
      #mydf <- data.frame(bug, ShannonDiversity, speciesRichness, caseControl, sex, wbo, ticLocation, hemsize_s_ml, bmi_CAT)
      mydf <- data.frame(bug, ShannonDiversity, speciesRichness, caseControl, sex, wbo, ticLocation, bmi_CAT)
      
      sampleSizeCase[index] <- length(myT[myT$caseControl == 1, i])
      sampleSizeControl[index] <- length(myT[myT$caseControl == 0, i])
      
      myLm <- lm( myT[ , i] ~  myT$wbo, na.action = na.exclude)
      myAnova <- anova(myLm)
      pValuesWBO[index] <- myAnova$"Pr(>F)"[1]
      rSquaredWBO[index] <- summary(myLm)$r.squared
   
      pValuesWBOKruskalWallis[index] <- pvalue(kruskal_test(bug ~ as.factor(wbo), data = mydf))
      # TODO Correct this 
      ## eta^2 is the effect size
      calcF <- qf(1 - pValuesWBOKruskalWallis[index], 3 - 1, length(bug) - 3)
      effectSizeWBOKruskalWallis[index] <- (calcF * (3 - 1))/(calcF * (3 - 1) + (length(bug) - 3))
      
      myLmBMICAT <- lm( myT[ , i] ~  myT$bmi_CAT, na.action = na.exclude)
      myAnovaBMICAT <- anova(myLmBMICAT)
      pValuesBMICAT[index] <- myAnovaBMICAT$"Pr(>F)"[1]
      rSquaredBMICAT[index] <- summary(myLm)$r.squared
      pValuesBMICATKruskalWallis[index] <- pvalue(kruskal_test(bug ~ as.factor(bmi_CAT), data = mydf))
      # TODO Correct this
      calcF <- qf(1 - pValuesBMICATKruskalWallis[index], 4 - 1, length(bug) - 4)
      effectSizeBMICATKruskalWallis[index] <- (calcF * (4 - 1))/(calcF * (4 - 1) + (length(bug) - 4))
      
      
      myLmWaist <- lm(myT[ , i] ~ myT$waist, na.action = na.exclude)
      myAnovaWaist <- anova(myLmWaist)
      pValuesWaist[index] <- myAnovaWaist$"Pr(>F)"[1]
      rSquaredWaist[index] <- summary(myLmWaist)$r.squared
      
      pValuesWaistKendall[index] <- Kendall(myT[ , i],  myT$waist)$sl[1]
      effectSizeWaistKendall[index] <- Kendall(myT[ , i],  myT$waist)$tau[1]
        #cor(myT[ , i], myT$waist, use = "na.or.complete", method = "kendall")
      
      myLmWHR <- lm(myT[ , i] ~ myT$whr, na.action = na.exclude)
      myAnovaWHR <- anova(myLmWHR)
      pValuesWHR[index] <- myAnovaWHR$"Pr(>F)"[1]
      rSquaredWHR[index] <- summary(myLmWHR)$r.squared
      
      pValuesWHRKendall[index] <- Kendall(myT[ , i],  myT$whr)$sl[1]
      effectSizeWHRKendall[index] <- Kendall(myT[ , i],  myT$whr)$tau[1]
      
      myLmAge <- lm(myT[ , i] ~ myT$age, na.action = na.exclude)
      myAnovaAge <- anova(myLmAge)
      pValuesAge[index] <- myAnovaAge$"Pr(>F)"[1]
      rSquaredAge[index] <- summary(myLmAge)$r.squared
      
      pValuesAgeKendall[index] <- Kendall(myT[ , i], myT$age)$sl[1]
      effectSizeAgeKendall[index] <- Kendall(myT[ , i], myT$age)$tau[1]
        #cor(myT[,i], myT$age, use="na.or.complete", method="kendall")
      #rSquaredAge[index] = rSquaredAge[index] * rSquaredAge[index]
      
      myLmBMI <- lm(myT[ , i] ~ myT$bmi, na.action = na.exclude)
      myAnovaBMI <- anova(myLmAge)
      pValuesBMI[index] <- myAnovaBMI$"Pr(>F)"[1]
      rSquaredBMI[index] <- summary(myLmBMI)$r.squared
      
      pValuesBMIKendall[index] <- Kendall(myT[ , i], myT$bmi)$sl[1]
      effectSizeBMIKendall[index] <- Kendall(myT[ , i], myT$bmi)$tau[1]
        #cor( myT[,i], myT$bmi, use="na.or.complete", method="kendall")
      #rSquaredBMI[index] = rSquaredBMI[index] * rSquaredBMI[index]

      #rSquaredValuesWaist[index] = rSquaredValuesWaist[index] * rSquaredValuesWaist[index]
      
      myLmTicsCount <- lm(myT[ , i] ~ myT$ticsCount, na.action = na.exclude)
      myAnovaTicsCount <- anova(myLmTicsCount)
      pValuesTicsCount[index] <- myAnovaTicsCount$"Pr(>F)"[1]
      rSquaredTicsCount[index] <- summary(myLmTicsCount)$r.squared
      
      pValuesTicsCountKendall[index] <- Kendall(myT[ , i] ,  myT$ticsCount)$sl[1] 
      effectSizeTicsCountKendall[index] <- Kendall(myT[ , i],  myT$ticsCount)$tau[1]
      #<- cor( myT[,i], myT$ticsCount, use="na.or.complete", method="kendall")
      #rSquaredTicsCount[index] = rSquaredTicsCount[index] * rSquaredTicsCount[index]
      
      #myLmCC <- lm(bug ~ as.factor(caseControl), na.action=na.exclude)
      #pValuesCaseControl[index] <- anova(myLmCC)$"Pr(>F)"[1]
      #rSquaredCaseControl[index] <- summary(myLmCC)$r.squared
      #pValuesCaseControl[index] <- 
      
      #  wilcox.test( myT[myT$caseControl==0,i] , myT[myT$caseControl==1,i])$p.value
      #pValuesCaseControl[index] <- wilcox_test(bug ~ caseControl, data = mydf, type="standardized")
      pValuesCaseControlWilcox[index] <- pvalue(wilcox_test(bug ~ as.factor(caseControl), data = mydf))
      rSquaredCaseControlWilcox[index] <- abs(statistic(wilcox_test(bug ~ as.factor(caseControl), data = mydf), type = "standardized")) / sqrt(length(bug))
      
      meanCase[index] <- mean(myT[myT$caseControl == 1, i])
      meanControl[index] <- mean(myT[myT$caseControl == 0, i])  
      
      myLm2 <- lm(myT[ , i] ~  factor(myT$caseControl), na.action = na.exclude)
      pValuesCaseControl[index] <- anova(myLm2)$"Pr(>F)"[1]
      rSquaredCaseControl[index] <- summary(myLm2)$r.squared
      
      #pValuesSex[index] <- 
      #  wilcox.test( myT[myT$sex == 1,i], myT[myT$sex == 2, i])$p.value
      
      pValuesSexWilcox[index] <- pvalue(wilcox_test(bug ~ as.factor(sex), data = mydf))
      rSquaredSexWilcox[index] <- abs(statistic(wilcox_test(bug ~ as.factor(sex), data = mydf), type="standardized"))/sqrt(length(bug))
      
      myLm3 <- lm( myT[,i] ~  factor(myT$sex) , na.action = na.exclude)
      pValuesSex[index] <- anova(myLm3)$"Pr(>F)"[1]
      rSquaredSex[index] <- summary(myLm3)$r.squared
      
      myLm4 <- lm( 
        myT[myT$ticLocation == "C__Distal_Only" | myT$ticLocation == "B__Proximal_Only", i] 
        ~  factor(myT[myT$ticLocation == "C__Distal_Only" | myT$ticLocation == "B__Proximal_Only", 11]) , 
        na.action=na.exclude)
      rSquaredticLocation[index] <- summary(myLm4)$r.squared
      pValuesticLocation[index] <- t.test(myT[myT$ticLocation == "C__Distal_Only", i],
                                         myT[myT$ticLocation == "B__Proximal_Only", i])$p.value
      
      storedf <- mydf
      mydf <- mydf[mydf$ticLocation == "C__Distal_Only" | mydf$ticLocation == "B__Proximal_Only",]
      
      #print(levels(mydf$ticLocation))
      mydf$ticLocation <- droplevels(mydf$ticLocation)
      pValuesTicLocationWilcox[index] <- pvalue(wilcox_test(mydf$bug ~ as.factor(mydf$ticLocation), data = mydf))
      rSquaredTicLocationWilcox[index] <- abs(statistic(wilcox_test(mydf$bug ~ as.factor(mydf$ticLocation), data = mydf), type="standardized"))/sqrt(length(mydf$bug))
      #print(dim(mydf))
      mydf <- storedf
      
      #myLm5 <- lm( 
      #  myT[myT$hemsize_s_ml=="Small" | myT$hemsize_s_ml=="Medium/Large",i] 
      #  ~  factor(myT[myT$hemsize_s_ml=="Small" | myT$hemsize_s_ml=="Medium/Large",11]) , 
      #  na.action=na.exclude)
      #rSquaredHemSize[index] <- summary(myLm5)$r.squared
      #pValuesHemSize[index] <- t.test(myT[myT$hemsize_s_ml=="Small",i],
      #                                myT[myT$hemsize_s_ml=="Medium/Large",i])$p.value
      
      #storedf <- mydf
      #mydf <- mydf[mydf$hemsize_s_ml != "none",]
      #mydf$hemsize_s_ml <- droplevels(mydf$hemsize_s_ml)
      #pValuesHemSizeWilcox[index] <- pvalue(wilcox_test(bug ~ as.factor(hemsize_s_ml), data = mydf))
      #rSquaredHemSizeWilcox[index] <- abs(statistic(wilcox_test(bug ~ as.factor(hemsize_s_ml), data = mydf), type="standardized"))/sqrt(length(bug))
      #print(dim(mydf))
      #mydf <- storedf
      
      names[index] <- names(myT)[i]
      
      boxplot(myT[myT$caseControl == 0, i], myT[myT$caseControl == 1, i], 
               main=paste(names[index], " p=", format(pValuesCaseControl[index], digits = 3) ))
      
      bug <- myT[,i]
      caseControl <- factor(myT$caseControl)
      myFrame <- data.frame(bug, caseControl)
      
      stripchart(bug ~ caseControl, 
                 data = myFrame, vertical = TRUE, pch = 21, add=TRUE, ylab = names[index])	
      #print(c("Makes it past the first boxplot", index))
      plot(myT$ticsCount, myT[ , i], main = paste(names[index], " p=", pValuesTicsCount[index]))
      plot(myT$waist, myT[ , i], main = paste(names[index], " p=", pValuesWaist[index]))
      plot(myT$whr, myT[ , i],  main = paste(names[index], " p=", pValuesWHR[index]))
      #print(c("Makes it past the first page of boxplots"))
      
      plot(myT[ , i] ~ myT$sex, main = paste(names[index], " p=", pValuesSex[index]))
      plot(myT[ , i] ~ myT$wbo, main = paste(names[index], " p=", pValuesWBO[index]))
      plot(myT[ , i] ~ myT$ticLocation, main = paste(names[index], " p=", pValuesticLocation[index]))
      #plot(myT[ , i] ~ myT$hemsize_s_ml, main = paste(names[index], " p=", pValuesHemSize[index]))
      plot.new()
      
      plot(myT$bmi, myT[ , i], main = paste(names[index], " p=", pValuesBMI[index]))
      plot(myT[ , i] ~ myT$bmi_CAT, main = paste(names[index], " p=", pValuesBMICAT[index]))
      plot(myT$age, myT[ , i], main = paste(names[index], " p=", pValuesAge[index]))
      plot.new()
      
      index <- index + 1
    }
  
  hist(pValuesCaseControl, breaks = 20)
  hist(pValuesCaseControlWilcox, breaks = 20)
  hist(pValuesSex, breaks = 20)
  hist(pValuesSexWilcox, breaks = 20)
  
  hist(pValuesticLocation, breaks = 20)
  hist(pValuesTicLocationWilcox, breaks = 20)
  hist(pValuesTicsCount, breaks = 20)
  hist(pValuesTicsCountKendall, breaks = 20)
  
  hist(pValuesWaist, breaks = 20)
  hist(pValuesWaistKendall, breaks = 20)
  hist(pValuesWHR, breaks = 20)
  hist(pValuesWHRKendall, breaks = 20)
  
  hist(pValuesWBO, breaks = 20)
  hist(pValuesWBOKruskalWallis, breaks = 20)
  hist(pValuesAge, breaks = 20)
  hist(pValuesAgeKendall, breaks = 20)
  
  hist(pValuesBMI, breaks=20)
  hist(pValuesBMIKendall, breaks = 20)
  hist(pValuesBMICAT, breaks = 20)
  hist(pValuesBMICATKruskalWallis, breaks = 20)
  
  #hist(pValuesHemSize, breaks = 20)
  #hist(pValuesHemSizeWilcox, breaks = 20)
  #plot.new()
  #plot.new()
  
  dev.off()
  
  pdf(paste("boxplots_Shannon_Richness_", taxa, ".pdf", sep = ""))
  par(mfrow = c(2, 2))
  

  
  myLm6 <- lm(ShannonDiversity ~ myT$caseControl)
  rSquaredShannon <- summary(myLm6)$r.squared
  pValuesShannon <- t.test(ShannonDiversity[myT$caseControl == 0],
                                  ShannonDiversity[myT$caseControl == 1])$p.value
  #pValuesWilcoxShannon <- wilcox.test( ShannonDiversity[myT$caseControl==0] , ShannonDiversity[myT$caseControl==1])$p.value
  pValuesWilcoxShannon <- pvalue(wilcox_test(ShannonDiversity ~ as.factor(caseControl), data = mydf))
  rSquaredWilcoxShannon <- abs(statistic(wilcox_test(ShannonDiversity ~ as.factor(caseControl), data = mydf), type="standardized")[1]) / sqrt(length(ShannonDiversity))
  #rSquaredShannonWilcox
  
  myLm7 <- lm(speciesRichness ~ myT$caseControl)
  rSquaredRichness <- summary(myLm7)$r.squared
  pValuesRichness <- t.test(speciesRichness[myT$caseControl == 0],
                           speciesRichness[myT$caseControl == 1])$p.value
  #pValuesWilcoxRichness <- wilcox.test( speciesRichness[myT$caseControl==0] , speciesRichness[myT$caseControl==1])$p.value
  pValuesWilcoxRichness <- pvalue(wilcox_test(speciesRichness ~ as.factor(caseControl), data = mydf))
  rSquaredWilcoxRichness <- abs(statistic(wilcox_test(speciesRichness ~ as.factor(caseControl), data = mydf), type="standardized")[1]) / sqrt(length(speciesRichness))
  #rSquaredRichnessWilcox
  
  ## sex
  
  ## bmi
  
  ## wbo
  
  ## bmi_CAT
  
  ## ticLocation
  
  
  ## hemsize_s_ml
  
  pValuesAgeRichness <- Kendall(speciesRichness, myT$age)$sl[1]
  rSquaredAgeRichness <- Kendall(speciesRichness, myT$age)$tau[1]
    #cor(speciesRichness, myT$age, use = "na.or.complete", method = "kendall")
  #rSquaredAgeRichness = rSquaredAgeRichness * rSquaredAgeRichness
  
  pValuesAgeShannon <- Kendall(ShannonDiversity, myT$age)$sl[1]
  rSquaredAgeShannon <- Kendall(speciesRichness, myT$age)$tau[1]
    #cor(ShannonDiversity, myT$age, use = "na.or.complete", method = "kendall")
  #rSquaredAgeShannon = rSquaredAgeShannon * rSquaredAgeShannon
  
  pValuesWaistRichness <- Kendall(speciesRichness, myT$waist)$sl[1]
  rSquaredWaistRichness <- Kendall(speciesRichness, myT$waist)$tau[1]
    #cor( speciesRichness, myT$waist, use="na.or.complete", method = "kendall")
  #rSquaredWaistRichness = rSquaredWaistRichness * rSquaredWaistRichness
  
  pValuesWaistShannon <- Kendall(ShannonDiversity, myT$waist)$sl[1]
  rSquaredWaistShannon <- Kendall(ShannonDiversity, myT$waist)$tau[1]
    #cor( speciesRichness, myT$waist, use="na.or.complete",method="kendall")
  #rSquaredWaistShannon = rSquaredWaistShannon * rSquaredWaistShannon
  
  pValuesWHRRichness <- Kendall(speciesRichness, myT$whr)$sl[1]
  rSquaredWHRRichness <- Kendall(speciesRichness, myT$whr)$tau[1]
  
  pValuesWHRShannon <- Kendall(ShannonDiversity, myT$whr)$sl[1]
  rSquaredWHRShannon <- Kendall(ShannonDiversity, myT$whr)$tau[1]
  
  pValuesTicsCountRichness <- anova(lm(speciesRichness ~ myT$ticsCount))$"Pr(>F)"[1]
  rSquaredTicsCountRichness <- summary(lm(speciesRichness ~ myT$ticsCount))$r.squared
  
  
  pValuesTicsCountRichnessKendall <- Kendall(speciesRichness, myT$ticsCount)$sl[1] 
  rSquaredTicsCountRichnessKendall <- Kendall(speciesRichness, myT$ticsCount)$tau[1]
    #cor( speciesRichness, myT$ticsCount, use="na.or.complete", method="kendall")
  #rSquaredTicsCountRichness = rSquaredTicsCountRichness * rSquaredTicsCountRichness
  
  pValuesTicsCountShannon <- anova(lm(ShannonDiversity ~ myT$ticsCount))$"Pr(>F)"[1]
  rSquaredTicsCountShannon <- summary(lm(ShannonDiversity ~ myT$ticsCount))$r.squared
  
  pValuesTicsCountShannonKendall <- Kendall(ShannonDiversity,  myT$ticsCount)$sl[1] 
  rSquaredTicsCountShannonKendall <- Kendall(ShannonDiversity, myT$ticsCount)$tau[1]
    #cor( ShannonDiversity, myT$ticsCount, use="na.or.complete",method="kendall")
  #rSquaredTicsCountShannon = rSquaredTicsCountShannon * rSquaredTicsCountShannon
  
  #caseControl = factor(caseControl, levels(caseControl)[c(2,1)])
  
  boxplot(ShannonDiversity ~ caseControl, 
          data = mydf, ylab="Shannon Diversity")
  
  stripchart(ShannonDiversity ~ caseControl, 
             data = mydf, vertical = TRUE, pch = 21, add=TRUE)	
  
  boxplot(speciesRichness ~ caseControl, 
          data = mydf)
  
  stripchart(speciesRichness ~ caseControl, 
             data = mydf, vertical = TRUE, pch = 21, add=TRUE)	
  
  plot(myT$ticsCount, ShannonDiversity)
  plot(myT$ticsCount, speciesRichness)
  
  dev.off()
  
  dFrame <- data.frame(names, sampleSizeCase, sampleSizeControl, meanCase, meanControl, 
                       pValuesCaseControl, rSquaredCaseControl, pValuesCaseControlWilcox, rSquaredCaseControlWilcox, 
                       pValuesSex, rSquaredSex, pValuesSexWilcox, rSquaredSexWilcox,
                       pValuesticLocation, rSquaredticLocation, pValuesTicLocationWilcox, rSquaredTicLocationWilcox,
                       #pValuesHemSize, rSquaredHemSize, pValuesHemSizeWilcox, rSquaredHemSizeWilcox,
                       pValuesWBO, rSquaredWBO, pValuesWBOKruskalWallis, effectSizeWBOKruskalWallis,
                       pValuesWaist, rSquaredWaist, pValuesWaistKendall, effectSizeWaistKendall,
                       pValuesWHR, rSquaredWHR, pValuesWHRKendall, effectSizeWHRKendall,
                       pValuesAge, rSquaredAge, pValuesAgeKendall, effectSizeAgeKendall,
                       pValuesBMI, rSquaredBMI, pValuesBMIKendall, effectSizeBMIKendall,
                       pValuesBMICAT, rSquaredBMICAT, pValuesBMICATKruskalWallis, effectSizeBMICATKruskalWallis, 
                       pValuesTicsCount, rSquaredTicsCount, pValuesTicsCountKendall, effectSizeTicsCountKendall)
  dFrameShannonRichness <- data.frame(pValuesShannon, rSquaredShannon, pValuesWilcoxShannon, rSquaredWilcoxShannon,
                       pValuesRichness, rSquaredRichness, pValuesWilcoxRichness, rSquaredWilcoxRichness,
                       pValuesTicsCountShannon, rSquaredTicsCountShannon, pValuesTicsCountShannonKendall, rSquaredTicsCountShannonKendall,
                       pValuesTicsCountRichness, rSquaredTicsCountRichness, pValuesTicsCountRichnessKendall, rSquaredTicsCountRichnessKendall,
                       pValuesWaistShannon, rSquaredWaistShannon, pValuesWaistRichness, rSquaredWaistRichness, 
                       pValuesAgeShannon, rSquaredAgeShannon, pValuesAgeRichness, rSquaredAgeRichness) 
  dFrame <- dFrame[order(dFrame$pValuesCaseControl),]
  dFrame$caseControlAdjust <- p.adjust(dFrame$pValuesCaseControl, method = "BH")
  dFrame$caseControlWilcoxAdjust <- p.adjust(dFrame$pValuesCaseControlWilcox, method = "BH")
  dFrame$pValuesSexAdjust <- p.adjust(dFrame$pValuesSex, method = "BH")
  dFrame$pValuesSexWilcoxAdjust <- p.adjust(dFrame$pValuesSexWilcox, method = "BH")
  # Might cause problems with subsetting
  dFrame$pValuesticLocationAdjust <- p.adjust(dFrame$pValuesticLocation, method = "BH")
  dFrame$pValuesTicLocationWilcoxAdjust <- p.adjust(dFrame$pValuesTicLocationWilcox, method = "BH")
  # Might cause problems with subsetting
  #dFrame$pValuesHemSizeAdjust <- p.adjust(dFrame$pValuesHemSize, method = "BH")
  #dFrame$pValuesHemSizeWilcoxAdjust <- p.adjust(dFrame$pValuesHemSizeWilcox, method = "BH")
  
  dFrame$pValuesWBOAdjust <- p.adjust(dFrame$pValuesWBO, method = "BH")
  dFrame$pValuesWBOKruskalWallisAdjust <- p.adjust(dFrame$pValuesWBOKruskalWallis, method = "BH")
  
  dFrame$pValuesWaistAdjust <- p.adjust(dFrame$pValuesWaist, method = "BH")
  dFrame$pValuesWaistKendallAdjust <- p.adjust(dFrame$pValuesWaistKendall, method = "BH")
  dFrame$pValuesWHRAdjust <- p.adjust(dFrame$pValuesWHR, method = "BH")
  dFrame$pValuesWHRKendallAdjust <- p.adjust(dFrame$pValuesWHRKendall, method = "BH")
  
  dFrame$pValuesTicsCountAdjust <- p.adjust(dFrame$pValuesTicsCount, method = "BH")
  dFrame$pValuesTicsCountKendallAdjust <- p.adjust(dFrame$pValuesTicsCountKendall, method = "BH")
  dFrame$pValuesAgeAdjust <- p.adjust(dFrame$pValuesAge, method = "BH")
  dFrame$pValuesAgeKendallAdjust <- p.adjust(dFrame$pValuesAgeKendall, method = "BH")
  dFrame$pValuesBMIAdjust <- p.adjust(dFrame$pValuesBMI, method = "BH")
  dFrame$pValuesBMIKendallAdjust <- p.adjust(dFrame$pValuesBMIKendall, method = "BH")
  dFrame$pValuesBMICATAdjust <- p.adjust(dFrame$pValuesBMICAT, method = "BH")
  dFrame$pValuesBMICATKruskalWallisAdjust <- p.adjust(dFrame$pValuesBMICATKruskalWallis, method = "BH")

  write.table(dFrame, file = paste("MCBT_metapValuesFor_", taxa, "_read1_.txt", sep = ""), sep = "\t", row.names = FALSE)
  write.table(dFrameShannonRichness, file = paste("MCBT_ShannonRichnessFor_", taxa, "_read1_.txt", sep = ""), sep = "\t", row.names = FALSE)
  
  # dFrameTicLocation <- data.frame(pValuesticLocation, rSquaredticLocation, pValuesTicLocationWilcox, rSquaredTicLocationWilcox)
  # dFrameTicLocation$pValuesticLocationAdjust <- p.adjust(dFrameTicLocation$pValuesticLocation, method = "BH")
  # dFrameTicLocation$pValuesticLocationWilcoxAdjust <- p.adjust(dFrameTicLocation$pValuesticLocationWilcox, method = "BH")
  # write.table(dFrameTicLocation, file = paste("MCBT_TicLocationpValuesFor_", taxa, "_read1_.txt", sep = ""), sep = "\t", row.names = FALSE)
  # 
  # dFrameHemSize <- data.frame(pValuesHemSize, rSquaredHemSize, pValuesHemSizeWilcox, rSquaredHemSizeWilcox)
  # dFrameHemSize$pValuesHemSizeAdjust <- p.adjust(dFrameHemSize$pValuesHemSize, method = "BH")
  # dFrameHemSize$pValuesHemSizeWilcoxAdjust <- p.adjust(dFrameHemSize$pValuesHemSizeWilcox, method = "BH")
  # write.table(dFrameHemSize, file = paste("MCBT_HemSizepValuesFor_", taxa, "_read1_.txt", sep = ""), sep = "\t", row.names = FALSE)
}
