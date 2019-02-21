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
	if (length(grep("^-source3=", myStr))>0)
	{
	source3 <- substring(myStr, nchar("-source3=")+1)
	source(source3) 
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
selected_cell<-rownames(raw_reads[raw_reads$TotalReads>1000000,])
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
#jpeg("/volumes/lab/users/minhu/Projects/2018-09-26-SCNA-kaile/output_1M/final_result/heatmap/heatmap.jpg", height=600, width=800)
#heatmap.3(as.matrix(mat), dendrogram="none", distfun = dist, hclustfun = function(x) hclust(x, method='ward.D2'),Colv=NA, 
#          ColSideColors=CHR, Rowv=TRUE,notecol="black",col=palettes,breaks=breaks, symm=F,symkey=F,symbreaks=T,trace="none",
#          cexRow=0.5, plot.row.partition=TRUE)
#dev.off()

tryCatch(expr = { library("pheatmap")}, 
         error = function(e) { 
           source("https://bioconductor.org/biocLite.R")
           biocLite("pheatmap")}, 
         finally = library("pheatmap"))


tryCatch(expr = { library("RColorBrewer")}, 
         error = function(e) { 
          install.packages("RColorBrewer")}, 
         finally = library("RColorBrewer"))

data_col_anno<-as.data.frame(cbind(chr))
colnames(mat)<-rownames(data_col_anno)

data_col_anno$chr<-gsub("2","n",data_col_anno$chr)
data_col_anno$chr<-gsub("1","2n",data_col_anno$chr)

ann_colors = list(chr = c("n"="#BEBEBE", "2n"= "#000000"))


my_hclust_cell <- hclust(dist(mat), method = "ward.D2")

# load package

tryCatch(expr = { library("dendextend")}, 
         error = function(e) { 
           source("https://bioconductor.org/biocLite.R")
           biocLite("dendextend")}, 
         finally = library("dendextend"))



#library(dendextend)
my_cell_col <- cutree(tree = as.dendrogram(my_hclust_cell), k = 3)
my_cell_col<-as.data.frame(my_cell_col)
my_cell_col$my_cell_col<-gsub("1","cluster 1",my_cell_col$my_cell_col)
my_cell_col$my_cell_col<-gsub("2","cluster 2",my_cell_col$my_cell_col)
my_cell_col$my_cell_col<-gsub("3","cluster 3",my_cell_col$my_cell_col)
colnames(my_cell_col) <- "clusters"

cv <- function(prof) sd(prof) / mean(prof)

# Input normalized ratio values
five.step.autocorrelation <- function(prof)
{
  right <- prof[-(1:5)]
  left <- prof[-((length(prof)-4):length(prof))]
  cor(left, right)
}


all<- data.frame()
for(x in 4:ncol(cnv_ratio)){
  
  temp<-cnv_ratio[,x]
  temp_cv<-cv(temp)
  all[(x-3),1]<-colnames(cnv_ratio)[x]
  all[(x-3),2]<-temp_cv
  temp_autoco<-five.step.autocorrelation(temp)
  all[(x-3),3]<-temp_autoco
  #all[(x-3),4]<-total_reads_per_cell2[all[(x-3),1],2]
  
}
rownames(all)<-all[,1]
rownames(all) <-gsub("\\.","-",rownames(all))
colnames(all)<-c("cell_id","cv","autoco")
all$cell_id<-NULL

all<-all[selected_cell,]
#all<-NaRV.omit(all)


my_cell_col$Total_Raw_Reads<-raw_reads[match(rownames(my_cell_col),rownames(raw_reads)),]$TotalReads
my_cell_col$Correlation_of_variation<-all[match(rownames(my_cell_col),rownames(all)),]$cv
my_cell_col$Autocorrelation_of_variation<-all[match(rownames(my_cell_col),rownames(all)),]$autoco
  
my_cell_col$Correlation_of_variation[which(my_cell_col$Correlation_of_variation>1)] =1
  

palettes <- colorRampPalette(c("blue", "white", "red"))(n = 999)
breaks = c(seq(-3,-0.5,length=50),seq(-0.49,-0.3,length=200),seq(-0.29,0.3,length=300),seq(0.31,0.5,length=400),seq(0.51,3,length=50)) 

heatmap_file<- paste(heatmap_dir,"heatmap_cluster.jpeg",sep="/")

jpeg(heatmap_file, height=1000, width=1200)
pheatmap(mat, annotation_row = my_cell_col, annotation_col = data_col_anno,annotation_colors=ann_colors,annotation_names_col = FALSE, color =palettes,breaks=breaks, clustering_method ="ward.D2",cluster_cols =FALSE,cutree_rows =3,show_colnames=FALSE)

dev.off()

heatmap_clusters <-paste(heatmap_dir,"heatmap.cluster.csv",sep="/")
write.csv(my_cell_col,file=heatmap_clusters)


