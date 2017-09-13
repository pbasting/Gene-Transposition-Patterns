#!/usr/bin/env Rscript

###############################################################################
#makeSyntenyPlot.r
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/1/2017
#Purpose: This is a component of a series of programs designed to classify protein
#		 'movement' when comparing two organisms and determine if proteins belonging
#		 to different functional categories are more likely to 'move'
#		 
#		 This function takes the results from 'CompareOrthologs' and makes a synteny
#		 plot with different colored rings indicating the different movement classifications
#		 	black = unmoved
#		 	blue = moved with adjacent proteins
#			red = moved without adjacent proteins
#		 	green = moved into a conserved region
#
#Arguments: (1)path to synteny results, (2)synteny result filename (no path), (3)name to give plot
################################################################################


args<-commandArgs(TRUE)
pathToFiles = args[1]
inputFile = args[2]
setwd(pathToFiles)
suppressMessages(require("lattice"))
suppressMessages(library(latticeExtra))

syntenyData <- read.csv(file=inputFile)


#parses out the data that corresponds to the proteins that moved and are no longer around the same proteins
get.Moved.Without.Adjacent <- function(syntenyData){
  subject <-c()
  query <-c()
  for(i in 1:nrow(syntenyData)){
    x = 1
    if (syntenyData[i, "Movement.Adjacent"] > 0){
      subject <- c(subject, syntenyData[i, "Subject.Protein"])
      query <- c(query, syntenyData[i, "Query.Protein"])
    }
  }
  moved.df <- data.frame(subject, query)
  return (moved.df)
}
movedWithoutAdjacent <- get.Moved.Without.Adjacent(syntenyData)

#parses out the data that corresponds to the proteins that moved alone and entered a conserved region
get.Moved.Without.Adjacent.Conserved <- function(syntenyData){
  subject <-c()
  query <-c()
  for(i in 1:nrow(syntenyData)){
    x = 1
    if (syntenyData[i, "Movement.Adjacent"] > 0){
       if (syntenyData[i, "Adjacent.Conserved"] > 0){
         subject <- c(subject, syntenyData[i, "Subject.Protein"])
         query <- c(query, syntenyData[i, "Query.Protein"])
      }
    }
  }
  moved.df <- data.frame(subject, query)
  return (moved.df)
}
conserved <- get.Moved.Without.Adjacent.Conserved(syntenyData)


#plot containing proteins that didn't move
completePlot <- xyplot(Subject.Protein ~ Query.Protein, data = syntenyData, col = "black", xlim=c(0,max(syntenyData$Query.Protein)),ylim=c(0,max(syntenyData$Subject.Protein)) ,  grid=TRUE)

#plot containing proteins that moved without adjacent proteins
if (length(movedWithoutAdjacent$subject) > 0){
	withoutPlot <- xyplot(subject ~ query, data = movedWithoutAdjacent,col = "blue")
}

#plot containing proteins that moved into a conserved region
if (length(conserved$subject) > 0){
  conservedPlot <- xyplot(subject ~ query, data = conserved,col = "red")
}


#checks that the plots exist as some may have no proteins
#Prevents failure due to missing plots
if (exists("conservedPlot") && exists("withoutPlot")){
	combinedplots <-completePlot+withoutPlot+conservedPlot
}else if (exists("withoutPlot")){
	combinedplots <- completePlot+withoutPlot
}else{
	combinedplots <- completePlot
}

#writes combined plto to pdf
pdf(args[3])
combinedplots
non <- dev.off()
