#!/usr/bin/env Rscript

###############################################################################
#makeSyntenyPlot.r
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/6/2017
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
moveCategories <-c(colnames(kegData)[2], colnames(kegData)[3], colnames(kegData)[4], colnames(kegData)[5], colnames(kegData)[6])
moveCatTotal <- c(sum(kegData$UNMOVED), sum(kegData$MOVED.ABS),sum(kegData$MOVED.ADJ), sum(kegData$MOVED.CONS),sum(kegData$MUTUAL.CONS))
moveTotal.df <- data.frame(moveCategories, moveCatTotal)

#total protein matches
totalProteins <- sum(moveCatTotal)

#makes dataframe containing the total matches for each keg category
kegCategories = kegData$FUNCTION
kegCatTotal <- rowSums(kegData[2:6])
kegTotal.df <- data.frame(kegCategories, kegCatTotal)


#calculates the expected average number of matches for each movement category and keg category
#stores in dataframe 'valuesExpected'
UN.exp <- c()
ABS.exp <- c()
ADJ.exp <- c()
CONS.exp <- c()
MUT.exp <- c()
for(i in 1:nrow(kegData)){
  UN.exp <- c(UN.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[1])/totalProteins))
  ABS.exp <- c(ABS.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[2])/totalProteins))
  ADJ.exp <- c(ADJ.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[3])/totalProteins))
  CONS.exp <- c(CONS.exp, (kegTotal.df$kegCatTotal[i]*(moveTotal.df$moveCatTotal[4])/totalProteins))
  MUT.exp <- c(MUT.exp, (kegTotal.df$kegCatTotal[i]*moveTotal.df$moveCatTotal[5]/totalProteins))
}
valuesExpected <- data.frame(kegCategories, UN.exp, ABS.exp, ADJ.exp, CONS.exp, MUT.exp)

#uses poisson distribution to calculate the likelihood that the given number would occur
#uses the lower tail if the actual count is lower than the average
#uses the upper tail if the acutal count is higher than the average
poisson <- function(actual, avg){
  if(actual > avg){
    pval = ppois(actual, lambda=avg, lower=FALSE)
  }else{
    pval = ppois(actual, lambda=avg, lower=TRUE)
  }
  return(pval)
}

#calculates the p-values and stores them in a dataframe 'poissonValues'
UN.pois <- c()
ABS.pois <- c()
ADJ.pois <- c()
CONS.pois <- c()
MUT.pois <- c()
for(i in 1:nrow(valuesExpected)){
  UN.pois <- c(UN.pois,poisson(kegData$UNMOVED[i],valuesExpected$UN.exp[i]))
  ABS.pois <- c(ABS.pois,poisson(kegData$MOVED.ABS[i],valuesExpected$ABS.exp[i]))
  ADJ.pois <- c(ADJ.pois,poisson(kegData$MOVED.ADJ[i],valuesExpected$ADJ.exp[i]))
  CONS.pois <- c(CONS.pois,poisson(kegData$MOVED.CONS[i],valuesExpected$CONS.exp[i]))
  MUT.pois <- c(MUT.pois,poisson(kegData$MUTUAL.CONS[i],valuesExpected$MUT.exp[i]))
}

poissonValues <- data.frame(kegCategories, kegData$UNMOVED,valuesExpected$UN.exp,UN.pois,
                            kegData$MOVED.ABS, valuesExpected$ABS.exp, ABS.pois,
                            kegData$MOVED.ADJ, valuesExpected$ADJ.exp, ADJ.pois,
                            kegData$MOVED.CONS, valuesExpected$CONS.exp, CONS.pois,
                            kegData$MUTUAL.CONS, valuesExpected$MUT.exp, MUT.pois)



#this generates a table where significant changes are assigned '+' or '-' depending
#on if the actual count is less or more than expected
#values with p values < 0.001 are assigned '++' or '--'
conservedIn1 <- c()
conservedIn2 <- c()
for(i in 1:nrow(poissonValues)){
  if(poissonValues$CONS.pois[i] < 0.05){
    if(poissonValues$CONS.pois[i] < 0.001){
      if(poissonValues$kegData.MOVED.CONS[i] > poissonValues$valuesExpected.CONS.exp[i]){
        conservedIn1 <-c(conservedIn1, "++")
      }else{
        conservedIn1 <-c(conservedIn1, "--")
      }
    }else{
      if(poissonValues$kegData.MOVED.CONS[i] > poissonValues$valuesExpected.CONS.exp[i]){
        conservedIn1 <-c(conservedIn1, "+")
      }else{
        conservedIn1 <-c(conservedIn1, "-")
      }
    }
    
  }else{
    conservedIn1 <-c(conservedIn1, " ")
  }
  
  if(poissonValues$MUT.pois[i] < 0.05){
    if(poissonValues$MUT.pois[i] < 0.001){
      if(poissonValues$kegData.MUTUAL.CONS[i] > poissonValues$valuesExpected.MUT.exp[i]){
        conservedIn2 <-c(conservedIn2, "++")
      }else{
        conservedIn2 <-c(conservedIn2, "--")
      }
    }else{
      if(poissonValues$kegData.MUTUAL.CONS[i] > poissonValues$valuesExpected.MUT.exp[i]){
        conservedIn2 <-c(conservedIn2, "+")
      }else{
        conservedIn2 <-c(conservedIn2, "-")
      }
    }
  }else{
    conservedIn2 <-c(conservedIn2, " ")
  }
  
}
resultsFigure <- data.frame(poissonValues$kegCategories, conservedIn1, conservedIn2)



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
