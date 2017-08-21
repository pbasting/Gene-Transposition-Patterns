#!/usr/bin/env Rscript
args<-commandArgs(TRUE)
pathToFiles = args[1]
inputFile = args[2]
setwd(pathToFiles)
require("lattice")
library(latticeExtra)

#inputFile <- "synteny_results/campylobacter_jejuni_NCTC_11168_and_campylobacter_coli_CVM_N29710/subject_campylobacter_jejuni_NCTC_11168_query_campylobacter_coli_CVM_N29710_MovementResults.csv"
syntenyData <- read.csv(file=inputFile)

#parses out the data that corresponds to the proteins that moved but are still near the same proteins
get.Moved.With.Adjacent <- function(syntenyData){
  subject <-c()
  query <-c()
  for(i in 1:nrow(syntenyData)){
    x = 1
    if (syntenyData[i, "Movement.Distance"] > 0){
      subject <- c(subject, syntenyData[i, "Subject.Protein"])
      query <- c(query, syntenyData[i, "Query.Protein"])
    }
  }
  moved.df <- data.frame(subject, query)
  return (moved.df)
}

movedWithAdjacent <- get.Moved.With.Adjacent(syntenyData)

#parses out the data that corresponds to the proteins that moved and are no longer around the same proteins
get.Moved.Without.Adjacent <- function(syntenyData){
  subject <-c()
  query <-c()
  for(i in 1:nrow(syntenyData)){
    x = 1
    if (syntenyData[i, "Movement.Distance"] > 0){
      if (syntenyData[i, "Movement.Adjacent"] > 0){
        subject <- c(subject, syntenyData[i, "Subject.Protein"])
        query <- c(query, syntenyData[i, "Query.Protein"])
      }
    }
  }
  moved.df <- data.frame(subject, query)
  return (moved.df)
}

movedWithoutAdjacent <- get.Moved.Without.Adjacent(syntenyData)

get.Moved.Without.Adjacent.Conserved <- function(syntenyData){
  subject <-c()
  query <-c()
  for(i in 1:nrow(syntenyData)){
    x = 1
    if (syntenyData[i, "Movement.Distance"] > 0){
      if (syntenyData[i, "Movement.Adjacent"] > 0){
        if (syntenyData[i, "Adjacent.Conserved"] > 0){
          subject <- c(subject, syntenyData[i, "Subject.Protein"])
          query <- c(query, syntenyData[i, "Query.Protein"])
        }
      }
    }
  }
  moved.df <- data.frame(subject, query)
  return (moved.df)
}

conserved <- get.Moved.Without.Adjacent.Conserved(syntenyData)

completePlot <- xyplot(Subject.Protein ~ Query.Protein, data = syntenyData, col = "black")

withPlot <- xyplot(subject ~ query, data = movedWithAdjacent, col = "blue")

withoutPlot <- xyplot(subject ~ query, data = movedWithoutAdjacent,col = "red")

conservedPlot <- xyplot(subject ~ query, data = conserved,col = "green")

combinedplots <-completePlot+withPlot+withoutPlot+conservedPlot

pdf(args[3])
combinedplots
dev.off()
