#!/usr/bin/env Rscript

#Segmentation of Profiles
#Author: Anna Unruh
#Last Updated: August 13, 2014

set.seed(1239)
theArgs <- NULL
for( myStr in commandArgs() )
{
  message("processing command arg:", myStr)
  if (length(grep("^-sample=", myStr))>0)
  {
    dataName <- substring(myStr, nchar("-sample=")+1) 
  }
  if (length(grep("^-undoprune=", myStr))>0)
  {
    undo.prune <- as.numeric(substring(myStr, nchar("-undoprune=")+1) )
  }
  if (length(grep("^-alpha=", myStr))>0)
  {
    alpha <- as.numeric(substring(myStr, nchar("-alpha=")+1))
  }
  if (length(grep("^-chrominfo=", myStr))>0)
  {
    chrominfo <- substring(myStr, nchar("-chrominfo=")+1)
  }
  if (length(grep("^-makeFig=", myStr))>0)
  {
    makeFig <- substring(myStr, nchar("-makeFig=")+1) 
  }
  if (length(grep("^-source1=", myStr))>0)
  {
    source1 <- substring(myStr, nchar("-source1=")+1) 
    source(source1)
  }
  if (length(grep("^-source2=", myStr))>0)
  {
    source2 <- substring(myStr, nchar("-source2=")+1) 
    source(source2)
  }
  if(length(grep("^-outdir=",myStr))>0)
  {
outArg <- substring(myStr, nchar("-outdir=")+1)
  }
}

print (source1)
print (source2)

dataDir <- paste(outArg, "blvb", sep="/")
outputDir <- paste(outArg, "final_result", sep="/")

print(dataName)
print(undo.prune)
print(alpha)
print(makeFig)
print(outputDir)
print(dataDir)
print(chrominfo)
if (!file.exists(outputDir)) {
  dir.create(outputDir)
}

CNProcessingCentromeres(dataName= dataName, dataDir = dataDir, outputDir = outputDir, centName=chrominfo,
                        alpha = alpha,undo.prune = undo.prune,
                        makeFig = makeFig)
