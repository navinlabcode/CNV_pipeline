CNProcessingCentromeres <- function(dataName = dataName,
                                    dataDir = dataDir,
                                    outputDir = outputDir,
                                    alpha = 0.01,
                                    min.width= 5,
                                    undo.prune = 0.05,
                                    fileName = ".vb",
                                    cpu = cpu,
                                    centName = centName,
                                    makeFig = FALSE)
                                    {
  set.seed(5471)

tryCatch(expr = { library("DNAcopy")},
         error = function(e) {
	if (!requireNamespace("BiocManager", quietly = TRUE))
    	install.packages("BiocManager")
        BiocManager::install("DNAcopy")},
        finally = library("DNAcopy"))
library(parallel)

  centr.loc <- read.delim(centName, as.is=TRUE)
  if (!file.exists(outputDir))
  { dir.create(outputDir) }

  if (!file.exists(paste(outputDir,"IGV_shorts", sep="/")))
  { dir.create(paste(outputDir,"IGV_shorts", sep="/")) }

  if (!file.exists(paste(outputDir,"IGV_longs", sep="/")))
  { dir.create(paste(outputDir,"IGV_longs", sep="/")) }

  cnBins <- NULL
  cnRatio <- NULL

  for(x in list.files(dataDir, pattern=fileName))
  {
    temp <- read.delim(paste(dataDir,x, sep="/"), header=FALSE)
    sample.name <- substring(x, 1, nchar(x) - nchar(fileName))
    sample.name <- sub(".bl","",sample.name)
    if(is.null(cnBins))
    {
      cnBins <- temp[,c(1:3,7)]
      cnRatio <- temp[,c(1:3,6)]
      colnames(cnBins)[ncol(cnBins)] <- sample.name
      colnames(cnRatio)[ncol(cnRatio)] <- sample.name
    } else
      if(!is.null(cnBins))
      {
        temp.bin <- temp$V7
        temp.r <- temp$V6
        cnBins <- cbind(cnBins, temp.bin)
        cnRatio <- cbind(cnRatio, temp.r)
        colnames(cnBins)[which(colnames(cnBins) == "temp.bin")] <- sample.name
        colnames(cnRatio)[which(colnames(cnRatio) == "temp.r")] <- sample.name
      }
  }

  colnames(cnBins)[match(c("V1", "V2", "V3"), colnames(cnBins))] <- c("chrom","chrompos","abspos")
  colnames(cnRatio)[match(c("V1", "V2", "V3"), colnames(cnRatio))] <- c("chrom","chrompos","abspos")

  write.table(cnBins, file = paste(outputDir, "/uber.", dataName, ".bin.txt", sep=""),
              quote = FALSE, sep = "\t",row.names = FALSE)
  write.table(cnRatio, file = paste(outputDir, "/uber.", dataName, ".ratio.txt", sep=""),
              quote = FALSE, sep = "\t",row.names = FALSE)
  
  cv <- function(n) sd(n) / mean(n)
  
  fun_parallel_seg<-function(x){
    sample.name <- colnames(cnRatio[,c(match(c("chrom","chrompos","abspos"), colnames(cnRatio)), x)])[4]
    thisRatio <-  cbind(cnBins[,c(match(c("chrom","chrompos","abspos"), colnames(cnBins)), x)],cnRatio[,x] )
    colnames(thisRatio) <-c("chrom", "chrompos", "abspos", "bincount", "ratio")
    thisRatio$ratio[which(thisRatio$ratio == 0)] <- 1e-3
    thisRatio$chrom <- as.numeric(thisRatio$chrom)
    CNA.object <- CNA(log(thisRatio$ratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio",
                      sampleid=sample.name)
    smoothed.CNA.object <- smooth.CNA(CNA.object)
        segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, min.width=5, undo.splits="prune", undo.prune=undo.prune)
    thisShort <- segment.smoothed.CNA.object[[2]]
    log.seg.mean.LOWESS <- rep(thisShort$seg.mean, thisShort$num.mark)
    merge.obj <- MergeLevels(log(thisRatio$ratio, base=2),log.seg.mean.LOWESS)

#keep absolution ratio data for plots
    seg.mean.LOWESS <- rep(2^thisShort$seg.mean, thisShort$num.mark)

# after mergeLevels
   log.seg.mean <- merge.obj$vecMerged
   seg.merge<- 2^merge.obj$vecMerged


################this step should be flexable for different datasets, default not to include this unless the samples are noisy

#  mean.log.seg <- median(log.seg.mean)
 # experimental.diff <- abs(log.seg.mean-log(thisRatio$ratio, base=2))
  # experimental.MAD <-median(experimental.diff)
   # m <- 1
    #for (m in 1:length(log.seg.mean)){
    #if(abs(log.seg.mean[m]-mean.log.seg)< 1.25*experimental.MAD){
   #log.seg.mean[m] <- mean.log.seg}
#m <- m+1
#}


##############################end of cutoff steps
meanseg <- 2^log.seg.mean
ratiox <- meanseg/mean(meanseg)
ratiox[ratiox==0] <- 0.0001
thisRatio$seg.mean.LOWESS <- ratiox

#making the short file after mergeLevels, note thisShort is written to seg file that keep raw segmented ratios
#thisRatio <- read.delim("/volumes/lab/ruligao/CNV_BRAC/06_13_2016_ES-BRCA-A_PM320/output_regular_parameter/IGV_longs/BRCA-A01.hg19.50k.varbin.txt")
#head(thisRatio)
#sample.name <- "test"

Short <- NULL
chr <- rle(thisRatio$chrom)[[2]]

for (x in 1:length(chr)){

	thisRatio.sub <- thisRatio[which(thisRatio$chrom==chr[x]), ]
	seg.mean.sub <- log(rle(thisRatio.sub$seg.mean.LOWESS)[[2]], base=2)
	rle.length.sub <- rle(thisRatio.sub$seg.mean.LOWESS)[[1]]

		num.mark.sub <- seq(1,length(rle.length.sub),1)
		loc.start.sub <-seq(1,length(rle.length.sub),1)
		loc.end.sub <- seq(1,length(rle.length.sub),1)

		len <-0
		j <-1

		for (j in 1: length(rle.length.sub)){
		num.mark.sub[j] <- rle.length.sub[j]
		loc.start.sub[j] <- thisRatio.sub$chrompos[len+1]
		len <- num.mark.sub[j]+len
		loc.end.sub[j] <- thisRatio.sub$chrompos[len]
		j <- j+1
		}

		ID <- rep(sample.name, times=length(rle.length.sub))
		chrom <- rep(chr[x], times=length(rle.length.sub))
		Short.sub <- cbind(ID,chrom,loc.start.sub,loc.end.sub,num.mark.sub,seg.mean.sub)
		Short <- rbind(Short, Short.sub)

		x <- x+1
}

colnames(Short) <- c("ID","chrom","loc.start","loc.end","num.mark","normalized_log2_seg.mean")

    if(makeFig)
    	{
 		 if (!file.exists(paste(outputDir,"Figures", sep="/")))
 		 { dir.create(paste(outputDir,"Figures", sep="/")) }

    	chr <- thisRatio$chrom
    	chr.shift <- c(chr[-1], chr[length(chr)])
    	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
    	hlines <- c(0.5, 1.0, 1.5, 2.0)
    	chr.text <- c(1:22, "X", "Y")
    	vlines.shift <- c(vlines[-1], 4*10^9)
    	chr.at <- vlines + (vlines.shift - vlines) / 2
    	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
    	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
    	y.at <- c(0.005, 0.02, 0.10, 0.50, 2.0, 6.00)
    	y.labels <- c("0.005", "0.02", "0.10", "0.50", "0.20", "6.00")
    	 print("Making figure")
 		jpeg(paste(outputDir,"/Figures/profile_", sample.name, ".jpg", sep=""), height=800, width=1200)

    	par(mar=c(5.1,4.1,6.1,4.1))
    	plot(x=thisRatio$abspos, y=thisRatio$ratio, log="y",ylim=range(0.05,6),main=c(paste(sample.name, ": GCnormalized+MergeLevels", sep=""),paste0("Average Read Counts Each Bin:  ",round(mean(thisRatio$bincount),0)),paste0("Coefficient of Variation (CV): ",round(cv(thisRatio$bincount),2))),xaxt="n", xlab="Genome Position (Gb)", yaxt="n", ylab="Copy Number Ratio", col="grey", cex=0.5)

    	axis(1, at=x.at, labels=x.labels)
    	axis(2, at=y.at, labels=y.labels)
    	lines(x=thisRatio$abspos, y=thisRatio$ratio, col="grey", cex=0.5)
        points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="blue",lwd=1.5, cex=0.5)
        lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="blue", lwd=1.5,cex=0.5)
    	#points(x=thisRatio$abspos, y=seg.mean.LOWESS, col="blue",lwd=1.5, cex=0.5)
    	#lines(x=thisRatio$abspos, y=seg.merge, col="red",lwd=1.5, cex=0.5)
        #points(x=thisRatio$abspos, y=seg.merge, col="red",lwd=1.5, cex=0.5)
        #lines(x=thisRatio$abspos, y=seg.mean.LOWESS, col="blue",lwd=1.5, cex=0.5)
    	abline(h=hlines, col="red")
    	abline(v= centr.loc$abs, col="red", lty = 2)
    	abline(v=vlines, par(lty=2))
    	mtext(chr.text, at = chr.at)
    	dev.off()

    }

write.table(thisRatio, sep="\t", file=paste(outputDir, "/IGV_longs/", sample.name, ".hg19.50k.varbin.txt",
                                            sep=""), quote=FALSE, row.names=FALSE)
write.table(Short, sep="\t", file=paste(outputDir, "/IGV_shorts/", sample.name, ".hg19.50k.varbin.txt",
                                        sep=""),  quote=FALSE, row.names=FALSE)

normlized_seg<-list(thisRatio$seg.mean.LOWESS)
names(normlized_seg)<-sample.name
sub_thisShort<-thisShort
names(sub_thisShort)<-c("ID","chrom","loc.start","loc.end","num.mark","Raw_log2_seg.mean")     
#return(normlized_seg)
combined_nor_log_seg_Short<-as.data.frame(Short)
names(combined_nor_log_seg_Short)<-c("ID","chrom","loc.start","loc.end","num.mark","Normalized_log2_seg.mean")   

results <- list(
  normlized_seg=normlized_seg,                             
  combined_log_thisShort=sub_thisShort,                             
  combined_nor_log_seg_Short=combined_nor_log_seg_Short   
)
return(results)
}
  
mc <- getOption("mc.cores", cpu)
res<-mclapply(4:ncol(cnRatio),fun_parallel_seg,mc.cores=mc)

all_seg<-as.data.frame(sapply(res,function(x){x$normlized_seg}))
chr_info<-cnRatio[, c("chrom","chrompos","abspos")]
all_seg<-cbind(chr_info,all_seg)

all_log_seg<-lapply(res,function(x){x$combined_log_thisShort})
all_log_seg<-do.call(rbind,all_log_seg)

all_nor_log_seg<-lapply(res,function(x){x$combined_nor_log_seg_Short})
all_nor_log_seg<-do.call(rbind,all_nor_log_seg)

write.table(all_seg, file = paste(outputDir, "/uber.", dataName, ".seg.txt", sep=""),
            quote = FALSE, sep = "\t",row.names = FALSE)
#  colnames(all_thisShort) <- c("ID","chrom","loc.start","loc.end","num.mark","Raw_seg.mean")

  write.table(all_log_seg, sep="\t", file=paste(outputDir, "/IGV_shorts/", dataName, ".all.hg19.50k.varbin.seg",
                                              sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)
  write.table(all_nor_log_seg, sep="\t", file=paste(outputDir, "/IGV_shorts/", dataName, ".all.hg19.50k.varbin.normalized.seg",
                                                sep=""), quote=FALSE, row.names=FALSE, col.names=TRUE)
}
