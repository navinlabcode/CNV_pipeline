facs.density <- function(inpath)
{
  # LOADING IN THE DATA
  # Load the FCS file as a flow frame
  flow.frame <- read.FCS(inpath)
  # Extract a tibble 
  flow.tibble.raw <- as_tibble(flow.frame@exprs)
  
  # From this tibble, extract a single numeric vector
  if (flow.frame@description$`$CYT`=="MoFlo Astrios")
  {
    raw.measurements <- flow.tibble.raw$`FL8-Area`
  } else if (flow.frame@description$`$CYT` == "FACSAriaII")
  {
    raw.measurements <- flow.tibble.raw$`DAPI-A`
  } else if (flow.frame@description$`$CYT` %in% c("BD FACSMelody", "FACSMelody"))
  {
    raw.measurements <- flow.tibble.raw$`DAPI*-A`
  } else
  {
    stop(sprintf('unrecognized flow cytometer "%s"',
                 flow.frame@description$`$CYT`))
  }
  # Remove stuff outside of range, by removing all maximum and minimum
  # observations. If all values are within range, two real observations will be
  # lost, which is acceptable
  measurements <- raw.measurements[raw.measurements!=max(raw.measurements) &
                                     raw.measurements!=min(raw.measurements)]
  
  # MAKING A KERNEL DENSITY ESTIMATE
  # Define the x values for the smoothed density
  xvals <- seq(min(raw.measurements), max(raw.measurements), length.out = 2^10)   #1024 values
  # We need to filter out NA values. I don't know why they're there in the first
  # place.
  intable <- data.frame(x = measurements[!is.na(measurements)])
  quadratic.model <- locfit(~ lp(x,nn=0.2,deg=2), data=intable)
  kernel.density.estimate <- 
    tibble(`DAPI fluorescence` = xvals,
           density=predict(quadratic.model, xvals))
  
  # DETECTING PEAKS IN THE KERNEL DENSITY ESTIMATE
  # Use the "peaks" function to detect peaks
  peaks.called <- kernel.density.estimate %>%
    mutate(is_peak = peaks(density))
  # Rank the peaks 
  # First, make a smaller table just containing the peaks, and rank them
  peak.ranks <- peaks.called %>%
    dplyr::filter(is_peak) %>%
    mutate(peak_rank=rank(-density, ties.method="min"))
  # Then, join it to the original, larger table containing both peaks and
  # non-peak values
  peaks.ranked <- left_join(
    peaks.called, dplyr::select(peak.ranks, -density, -is_peak),
    by="DAPI fluorescence"
  )
  return(peaks.ranked)
}

