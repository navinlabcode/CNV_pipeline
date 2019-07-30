#!/usr/bin/env Rscript

# Heatmap for all CNV Profiles
# Author: Min Hu
# Last Updated: Nov 3th, 2018

set.seed(100)

for (myStr in commandArgs() )
{
	message(" processing command arg:",myStr)
	if(length(grep ("^-sample=", myStr)) >0)
	{
	dataName <- substring(myStr, nchar("-sample=")+1) 
	}
	if(length(grep("^-outdir=",myStr))>0)
	{
	outdir <- substring(myStr, nchar("-outdir=")+1)
	}
}
#data
final_dir <- paste(outdir,"final_result",sep="/")
heatmap_dir<-paste(final_dir,"heatmap",sep="/")
if(!file.exists(heatmap_dir)){
 dir.create(heatmap_dir)
}

segfile=list.files(final_dir,pattern="seg.txt")

segmentInfo  <- read.delim(paste(final_dir,segfile,sep="/"), header= TRUE)
##uber input uber.seg.txt, which contains segmented ratio values in columns of each cell, first three columns contain chr, start, absolute position

##reads files
stat_dir<-paste(outdir,"metrics",sep="/")
stat_file<-list.files(stat_dir,pattern="all_stat_metrics.txt")
raw_reads<-read.delim(paste(stat_dir,stat_file,sep="/"),header=TRUE,sep="\t")
rownames(raw_reads)<-raw_reads$Sample.Name
raw_reads$Sample.Name<-NULL

##CNV ratio files

ratio_file=list.files(final_dir,pattern="ratio.txt")
cnv_ratio<-read.table(paste(final_dir,ratio_file,sep="/"),header=TRUE)


## load library
#library("devtools")
n<-ncol(segmentInfo)

##get the segment ratio data
sam <- t(as.matrix(segmentInfo[,4:n]))

##assign column names
colnames(sam) <- as.vector(segmentInfo[,1])

##asign sample names as rownames
rownames(sam) <- colnames(segmentInfo)[4:n]
rownames(sam) <-gsub("\\.","-",rownames(sam))

## filter out the cells with Reads counts less than 1M
selected_cell<-rownames(raw_reads[raw_reads$TotalReads>100000,])
sam<-sam[match(selected_cell,rownames(sam)),]
dim(sam)

tryCatch(expr = { library("flowViz")}, 
         error = function(e) { 
           source("https://bioconductor.org/biocLite.R")
           biocLite("flowViz")}, 
         finally = library("flowViz"))

#require(IDPmisc)
#sam<-NaRV.omit(sam)

# transform seg ratio data for analysis
mat <-as.matrix(sam) +0.1
mat <- log2(mat)


# hierarchycal clustering ############################################


#library(gplots)

###set colobreaks for heatmap
palettes <- colorRampPalette(c("blue", "white", "red"))(n = 999)
breaks = c(seq(-2,-0.5,length=50),seq(-0.5,-0.3,length=200),seq(-0.3,0.3,length=300),seq(0.3,0.5,length=400),seq(0.5,3,length=50)) 

##color bars
chr<- as.numeric(colnames(mat)) %% 2+1
rbpal <- colorRampPalette(c("grey","black"))
CHR <- cbind(rbpal(2)[as.numeric(chr)], rbpal(2)[as.numeric(chr)])



#####heatmap
heatmap_file<- paste(heatmap_dir,"heatmap_cluster.jpeg",sep="/")

jpeg(heatmap_file, height=1000, width=1200)
source("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
heatmap.3(as.matrix(mat), dendrogram="none", distfun = dist, hclustfun = function(x) hclust(x, method='ward.D2'),Colv=NA, 
          ColSideColors=CHR, Rowv=TRUE,notecol="black",col=palettes,breaks=breaks, symm=F,symkey=F,symbreaks=T,trace="none",
          cexRow=0.5, plot.row.partition=TRUE)
dev.off()

