library(tidyverse)
library(flsa)

## FUNCTIONS USED FOR CONVERTING BINCOUNTS INTO VALUES AND WEIGHTS

# Convert bincounts to logratios and segment
lrp <- function(x, extra=.15)
{
  log2(x / mean(x) + extra)
}
segmented.ratio.long <- function(corrected.bincounts, penalty=sqrt(10), extra=.15)
{
  flsaGetSolution(flsa(lrp(corrected.bincounts[!is.na(corrected.bincounts)], extra=extra)), lambda2=penalty)[1,]
}
# Functions for breaking up longs, summarizing as medians, calculating weights,
# and filtering out small segments
segnum <- function(long) cumsum(c(TRUE, diff(long)!=0))
# I know segsummary is just tapply but it might not be in the future, since you
# can do about ten times better than tapply using a foreign function call
segsummary <- function(unseg, these.segnums, sumfun) tapply(unseg, these.segnums, sumfun)
seglengths <- function(these.segnums) rle(these.segnums)$lengths
h <- 2^-6
summarize.seg <- function(gc_bincount, long_ratios, sumfun=median, required.size = 40)
{
  these.segnums <- segnum(long_ratios);
  unseg.ratios <- gc_bincount/mean(gc_bincount);
  summarized.prof <- tibble(ratios=segsummary(unseg.ratios, these.segnums, sumfun),
         lengths=seglengths(these.segnums))
  weighted.prof <- summarized.prof %>%
    mutate(weights=2*(pnorm(h, sd=sqrt(lengths/ratios))-1/2))
  filtered.prof <- weighted.prof %>%
    dplyr::filter(lengths >= required.size)
  return(filtered.prof)
}
# One function that does all of the above
summarize.bincounts <- function(gc.bincount, extra=.15, sumfun=median, penalty=sqrt(10), required.size = 40)
{
  long.ratios <-
    segmented.ratio.long(gc.bincount, penalty=penalty, extra=extra)
  sample.for.ploidy.inference <-
    summarize.seg(gc.bincount, long.ratios, sumfun=sumfun, required.size=required.size)
  return(sample.for.ploidy.inference)
}

## FUNCTIONS FOR PLOIDY INFERENCE
# Function to make the histogram
weighted.hist <- function(xs, ws)
{
  breaks <- c(seq(0,2^6, length.out=2^12), Inf)
  histogram <- tapply(ws, cut(xs, breaks), sum)
  histogram[is.na(histogram)] <- 0
  xvals <- (1:length(histogram))/2^6
  return(tibble(x=xvals,  y=histogram))
}
ecf <- function(histogram)
{
  # x values from the histogram are never used
  transformed <- fft(histogram$y, inverse=TRUE)/sum(histogram$y)
  frequencies <- (1:2^12)/(2^6)
  # Using tibble guarantees x and y will be the same length, or there will be an
  # error, or one of the values will be repeated
  return(tibble(x=frequencies, y=transformed))
}
weighted.ploidy.estimate <- function(values, lengths, weights)
{
  stopifnot(length(lengths)==length(values))
  stopifnot(length(lengths)==length(weights))
  if (length(lengths)==0) return(tibble(ploidy = NA,
         height = NA,
         second_height=NA))
  filtered.values <- values[!is.infinite(weights)]
  # Calculate ratios
  meanval <- sum(lengths*values) / sum(lengths)
  xs <- filtered.values / meanval
  ws <- weights[!is.infinite(weights)]
  histogram <- weighted.hist(xs, ws)
  stopifnot(nrow(histogram) == 2^12)
  this.ecf <- ecf(histogram)
  output <- with(this.ecf,
  {
    target.range <- x >= 1 & x <= 8
    # This assumes x and y are the same length, so a logical vector taken from x
    # can be used to index y
    # fx means "filtered x", likewise with fy
    fx <- x[target.range]
    fy <- y[target.range]
    top.index <- which.max(Mod(fy))
    # fx and fy were indexed by the same logical vector, and guaranteed to be
    # the same length. Therefore an index taken from fy can be used to pull a
    # value from fx
    tibble(ploidy = fx[top.index],
         height = max(Mod(fy)),
         second_height=max(Mod(fy[round(top.index*1.5):length(fy)])))
  })
  return(output)
}

## FUNCTION TO DO THE WHOLE PROCESS
ploidy.and.peakheight <- function(bincounts, extra=.15, sumfun=median, penalty=sqrt(10), required.size = 40)
{
  sample.for.ploidy.inference <- 
    summarize.bincounts(bincounts, extra=extra, sumfun=sumfun, penalty=penalty, required.size = required.size)
  results <- with(sample.for.ploidy.inference, weighted.ploidy.estimate(ratios, lengths, weights))
  return(with(results, c(ploidy=ploidy, peak_height=height)))
}