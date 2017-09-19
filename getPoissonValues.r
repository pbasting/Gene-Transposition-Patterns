#!/usr/bin/env Rscript

###############################################################################
#makeSyntenyPlot.r
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/19/2017
#Purpose: This is a component of a series of programs designed to classify protein
#		 'movement' when comparing two organisms and determine if proteins belonging
#		 to different functional categories are more likely to 'move'
#		 
#		 This function takes the results from 'FormatKegResults' and calculates the 
#		 average expected matches for a given keg category. A Poisson approximation
#		 is used to estimate the probability that a given number of matches would
#		 excede or be lower than the expected average. A table is generated to
#		 easily show which categories have staistically significant changes
#
#Arguments: (1)csv output from FormatKegResults (2)name of poisson results csv (3)name of summary table csv
################################################################################

args<-commandArgs(TRUE)
inputFile = args[1]
poisOutput = args[2]
figureOutput = args[3]

kegData <- read.csv(file=inputFile)

#makes dataframe containing total matches for each movement category
moveCategories <-c(colnames(kegData)[2], colnames(kegData)[3], colnames(kegData)[4], colnames(kegData)[5])
moveCatTotal <- c(sum(kegData$UNMOVED), sum(kegData$MOVED), sum(kegData$MOVED.CONS),sum(kegData$MUTUAL.CONS))
moveTotal.df <- data.frame(moveCategories, moveCatTotal)

#total protein matches
totalProteins <- sum(moveCatTotal)

#makes dataframe containing the total matches for each keg category
kegCategories = kegData$FUNCTION
kegCatTotal <- rowSums(kegData[2:5])
kegTotal.df <- data.frame(kegCategories, kegCatTotal)


#calculates the expected average number of matches for each movement category and keg category
#stores in dataframe 'valuesExpected'
UN.exp <- c()
MOV.exp <- c()
CONS.exp <- c()
MUT.exp <- c()
for(i in 1:nrow(kegData)){
  UN.exp <- c(UN.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[1])/totalProteins))
  MOV.exp <- c(MOV.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[2])/totalProteins))
  CONS.exp <- c(CONS.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[3])/totalProteins))
  MUT.exp <- c(MUT.exp, (kegTotal.df$kegCatTotal[i]*moveTotal.df$moveCatTotal[4]/totalProteins))
}
valuesExpected <- data.frame(kegCategories, UN.exp, MOV.exp, CONS.exp, MUT.exp)

#uses poisson distribution to calculate the likelihood that the given number would occur
#uses the lower tail if the actual count is lower than the average
#uses the upper tail if the acutal count is higher than the average

poisson <- function(actual, avg){
  dif = abs(actual - avg)
  pval = (ppois(avg+dif, lambda=avg, lower=FALSE))+(ppois(avg-dif,lambda=avg, lower=TRUE))
  return(pval)
}

#calculates the p-values and stores them in a dataframe 'poissonValues'
UN.pois <- c()
MOV.pois <- c()
CONS.pois <- c()
MUT.pois <- c()
for(i in 1:nrow(valuesExpected)){
  UN.pois <- c(UN.pois,poisson(kegData$UNMOVED[i],valuesExpected$UN.exp[i]))
  MOV.pois <- c(MOV.pois,poisson(kegData$MOVED[i],valuesExpected$MOV.exp[i]))
  CONS.pois <- c(CONS.pois,poisson(kegData$MOVED.CONS[i],valuesExpected$CONS.exp[i]))
  MUT.pois <- c(MUT.pois,poisson(kegData$MUTUAL.CONS[i],valuesExpected$MUT.exp[i]))
}

#This adjusts pvals to account for fasle discovery rate
allpvals <- c(UN.pois, MOV.pois, CONS.pois, MUT.pois)
padjBH <- p.adjust(allpvals, method = "BH")
for(i in 1:length(UN.pois)){
  UN.pois[i] <- padjBH[i]
}
for(i in 1:length(MOV.pois)){
  MOV.pois[i] <- padjBH[i+length(UN.pois)]
}
for(i in 1:length(CONS.pois)){
  CONS.pois[i] <- padjBH[i+length(UN.pois)+length(MOV.pois)]
}
for(i in length(MUT.pois)){
  MUT.pois[i] <- padjBH[i+length(UN.pois)+length(MOV.pois)+length(CONS.pois)]
}

poissonValues <- data.frame(kegCategories, kegData$UNMOVED,valuesExpected$UN.exp,UN.pois,
                            kegData$MOVED, valuesExpected$MOV.exp, MOV.pois,
                            kegData$MOVED.CONS, valuesExpected$CONS.exp, CONS.pois,
                            kegData$MUTUAL.CONS, valuesExpected$MUT.exp, MUT.pois)


#this generates a table where significant changes are assigned '+' or '-' depending
#on if the actual count is less or more than expected
#values with p values < 0.001 are assigned '++' or '--'
determineEnrichment <- function(poisson, actual, expected, column){
  if(poisson < 0.05){
    if(poisson < 0.001){
      if (actual > expected){
        column <-c(column, "++")
      }else{
        column <-c(column, "--")
      }
    }else{
      if(actual > expected){
        column <-c(column, "+")
      }else{
        column <-c(column, "-")
      }
    }
  }else{
  	column <- c(column, " ")
  }
  return(column)
}

unmoved <- c()
unconserved <- c()
conservedIn1 <- c()
conservedIn2 <- c()
for(i in 1:nrow(poissonValues)){
  unmoved <- determineEnrichment(poissonValues$UN.pois[i], poissonValues$kegData.UNMOVED[i], poissonValues$valuesExpected.UN.exp[i], unmoved)
  unconserved <- determineEnrichment(poissonValues$MOV.pois[i], poissonValues$kegData.MOVED[i], poissonValues$valuesExpected.MOV.exp[i], unconserved)
  conservedIn1 <- determineEnrichment(poissonValues$CONS.pois[i], poissonValues$kegData.MOVED.CONS[i], poissonValues$valuesExpected.CONS.exp[i], conservedIn1)
  conservedIn2 <- determineEnrichment(poissonValues$MUT.pois[i], poissonValues$kegData.MUTUAL.CONS[i], poissonValues$valuesExpected.MUT.exp[i], conservedIn2)
}
resultsFigure <- data.frame(poissonValues$kegCategories, unmoved, unconserved, conservedIn1, conservedIn2)



#finds the keg categories that have no matches
noMatches <- c()
for(i in 1:nrow(kegTotal.df)){
  if(kegTotal.df$kegCatTotal[i] == 0){
    noMatches <- c(noMatches, i)
  }
}

#removes the rows with no matches
poissonValues <- poissonValues[-noMatches,]
resultsFigure <- resultsFigure[-noMatches,]


#outputs results to files
write.csv(poissonValues, poisOutput)
write.csv(resultsFigure, figureOutput)
