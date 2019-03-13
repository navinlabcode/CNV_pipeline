#!/usr/bin/env Rscript

#Filtering vb files for Copy Number Data
#Description: Removes blacklisted files 
#Depends on: R.
#Requires: 
#Author: Anna Unruh
#Author: Min
#updated: correct the bugs with empty input file, zero-inflated and NA-inflated distributions
#Last Updated: Nov 13, 2018

theArgs <- NULL
for( myStr in commandArgs() )
{
	  message("processing command arg:", myStr)
  	if (length(grep("^-filter=", myStr))>0)
  	{
    	filterFile <- substring(myStr, nchar("-filter=")+1) 
  	}
  	if (length(grep("^-gcinputFile=", myStr))>0)
 	 {
    	gcInputFile <- substring(myStr, nchar("-gcinputFile=")+1) 
  	}
 	if(length(grep("^-outdir=",myStr))>0)
 	{
 	outArg <- substring(myStr, nchar("-outdir=")+1) 
 	}
 	if(length(grep("^-filter_CellWithEmptyBin=",myStr))>0)
 	{
 	filter_CellWithEmptyBin <-as.numeric(substring(myStr, nchar("-filter_CellWithEmptyBin=")+1))
	 }
	
}
#thisArgs <- read.delim(theArgs[1], header=FALSE, as.is=TRUE)
#thisArgs <- thisArgs$V1

#Read in the arguments from your arguement file`s
#filterFile <- thisArgs[1]
#varbinPattern <- thisArgs[2]
#gcInputFile <- thisArgs[3]
varbinPattern<-".vb"

vb_folder <- paste(outArg,"vbdir",sep = "/")
blvb <- paste(outArg,sub(".","bl",varbinPattern),sep = "/")

#Make sure that the files exist
if(!file.exists(filterFile))
	stop("The filter file you specified does not exist.")
print(vb_folder)	
print (varbinPattern)
#varbinPattern <- ".vb"
testFile <- list.files(vb_folder, pattern=varbinPattern, full.names=TRUE)[1]
print(testFile)
if(!file.exists(testFile))
	stop("The varbin directory shows no existing files with the extension you specified.")
if(!file.exists(gcInputFile))
	stop("The GC file you specified does not exist.")

#read in the files
filter <- read.delim(filterFile)
varb <- read.delim(list.files(vb_folder, pattern=varbinPattern, full.names=TRUE)[1], header=FALSE)
gcFile <- read.delim(gcInputFile, header=TRUE)

#Check that the first varbin file has the same number of rows as the GC file
if(nrow(varb) != nrow(gcFile))
	{
	stop(paste(sep="\n", "Varbins files do not have the same number of rows as GC files.", 
               paste("The varbin files have", nrow(varb), "rows. \nThe GC file has", nrow(gcFile), "rows."),
               "Check the varbin directory and files.",
               "The GC file you specified was:",
                gcInputFile))
        }
#LOWESS GC normalization
lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}


#Records a vector containing the chromomsome information for the varbin files
vChrNames <- sapply(varb$V1, function(x)
			{
			reValue <- x
			if(x == 23)
				{reValue <- "X"}
			if(x == 24)
				{reValue <- "Y"}
			reValue
				})

#checks that column names contain th start bin
if(!("bin.start" %in% colnames(filter)))
	stop("Filter file does not contain the correct column names.")
if(!("bin.chrom" %in% colnames(filter)))
	stop("Filter file does not contain the correct column names.")	

gc <- read.delim(gcInputFile, header=TRUE)

#Make sure we are removing 1 line for each of the to Remove lines (No more or less)
toRemove <- NULL
	if(grepl("2M|5M|10M",filterFile))
	{ 
		toRemove <- NULL
	gcFiltered <- gc
	}else{
	for(i in 1:nrow(filter))
		{
		if(length(which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & varb$V2 == (filter$bin.start[i] + 1))) == 1)
			{ 
			toRemove <- append(toRemove, 
				which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & 
				varb$V2 == (filter$bin.start[i] + 1)))
			} else 
			if(length(which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & varb$V2 == (filter$bin.start[i]))) == 1)
			{
			toRemove <- append(toRemove, 
				which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & 
				varb$V2 == (filter$bin.start[i])))		
			} else
			if(length(which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & varb$V2 != (filter$bin.start[i]))) == 1 | 
			length(which(paste("chr", vChrNames, sep="") == filter$bin.chrom[i]  & varb$V2 == (filter$bin.start[i] + 1))) == 1)
				{
				stop("Problem the columns in the bins in the varbin and gc files do not match.")
				}
		}
	
	gcFiltered <- gc[-toRemove,] #Make gc file - filtered
	}
#gc <- read.delim(gcInputFile, header=TRUE)
#gcFiltered <- gc[-toRemove,]


#Make varbins file filtered
print("Writing out the varbins files.")
if(!file.exists(blvb))
  { dir.create(blvb)}

for(file in list.files(vb_folder, pattern=varbinPattern, full.names=TRUE))
	{
      if(file.info(file)$size==0){
       	next
      }		
     vabtemp <- read.delim(file, header=FALSE)
	print (file)
	print (length(which(vabtemp$V5==0)))
	print (filter_CellWithEmptyBin*nrow(vabtemp))
 
     if(length(which(vabtemp$V5==0)) > filter_CellWithEmptyBin*nrow(vabtemp)){
      	next
     	 }	 
                vabtemp$V5[vabtemp$V5==0] <- 0.001
		vabtemp$V6 <- lowess.gc(gcFile$gc.content,vabtemp$V5)
                vabtemp$V7 <- vabtemp$V6*median(vabtemp$V4)
                if(sum(is.na(vabtemp$V7))>0){
		next	
		}
		temp <- strsplit(file, "/")[[1]]
		temp2 <- strsplit(temp[length(temp)],split=".", fixed=TRUE)[[1]]
		sampleName <- temp2[1]
	if(is.null(toRemove)){	
	write.table(x=vabtemp, file=paste(blvb,  paste(sampleName,"-bl", varbinPattern, sep=""), sep="/"), col.names=FALSE, row.names=FALSE, sep="\t")
	}else{
	write.table(x=vabtemp[-toRemove,], file=paste(blvb,  paste(sampleName,"-bl", varbinPattern, sep=""), sep="/"), col.names=FALSE, row.names=FALSE, sep="\t")
	}
}
