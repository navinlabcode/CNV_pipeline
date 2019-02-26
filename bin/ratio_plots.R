ratio_plots <- function(seg_data,
                        bin_data,
                        ncpu = 40) {
  
  # obtaining integer median value to use for color
  source("bin/integer_cn.R")
  
  #copying bin_data to calculate ploidy
  bin_data_cp <- bin_data
  if (any(names(bin_data_cp) %in% c("chrom",
                                 "chrompos",
                                 "abspos"))) {
    bin_data_cp <- dplyr::select(bin_data_cp, -chrom, -chrompos, -abspos) %>% 
      as.data.frame()
  }
  
  browser()
  
  # calling Alex integer method sourced from integer_cn.R
  dat_ploidy <- Reduce(rbind, 
                       parallel::mclapply(bin_data_cp, 
                                          ploidy.and.peakheight, 
                                          mc.cores= ncpu)) %>% 
    as_tibble()
  
  # saving the median
  ploidy_median <- median(dat_ploidy$ploidy)
  
  # establishing a value to truncate the ploidy that will be used for the colors later
  ploidy_trunc_val <- round(ploidy_median) * 2
  
    # Creating absolute position of chromosomes, code slighly modified from Alex's ggplot CNV heatmap
  chrom.lengths <-
    c(
      249250621,
      243199373,
      198022430,
      191154276,
      180915260,
      171115067,
      159138663,
      146364022,
      141213431,
      135534747,
      135006516,
      133851895,
      115169878,
      107349540,
      102531392,
      90354753,
      81195210,
      78077248,
      59128983,
      63025520,
      48129895,
      51304566,
      155270560,
      59373566
    )
  chrom.names <- c(paste("chr", 1:22, sep = ""),
                   "chrX",
                   "chrY")
  chrom.abspos.end <- cumsum(chrom.lengths)
  names(chrom.abspos.end) <- chrom.names
  chrom.abspos.start <-
    c(0, chrom.abspos.end[1:(length(chrom.abspos.end) - 1)])
  names(chrom.abspos.start) <- names(chrom.abspos.end)
  
  # Creating data frame object for chromosome rectangles shadows
  chrom.rects <- data.frame(xstart = chrom.abspos.start,
                            xend = chrom.abspos.end)
  xbreaks <- rowMeans(chrom.rects)
  chrom.rects$colors <-rep(c("white", "gray"), 12)
  
  #storing chromosomes info
  chrom_info <- bin_data %>% dplyr::select(c(abspos,chrom, chrompos))
  
  # removing chromosomes info to calculate the ratios
  bin_data <- bin_data %>% dplyr::select(-c(abspos,chrom, chrompos))
  ratios_data <- sweep(bin_data, 2, apply(bin_data, 2, mean), '/')
  
  # removing chromosomes info from seg_data
  seg_data <- seg_data %>% dplyr::select(-chrom, -chrompos)
  
  # adding chromosomes info back
  ratios_data$abspos <- chrom_info$abspos
  
  # ggplot will need the long format tables, using gather
  ratios_data_l_gat <- tidyr::gather(data = ratios_data, key = "cell", value = "ratios_rat", -abspos)
  long_seg_t_gat <- tidyr::gather(data = seg_data, key = "cell", value = "seg_mean", -abspos)
  # I'll need to join both tables
  df <- inner_join(ratios_data_l_gat, long_seg_t_gat)
  # ploidy scaling the seg values
  df <- df %>% dplyr::mutate(int_value = round(ploidy_median*seg_mean))
  # obtaining the max value. I'll use later on the plot
  max_int_value <- max(df$int_value)
  # truncating so I can limit the colors
  df$int_value[df$int_value > ploidy_trunc_val] <- ploidy_trunc_val
  
  # setting colors of the ratios according to the inferred ploidy
  if (round(ploidy_median) == 2) {
    df <- df %>% dplyr::mutate(colors = case_when(
      int_value == 0 ~ "darkblue",
      int_value == 1 ~ "dodgerblue",
      int_value == 2 ~ "azure4",
      int_value == 3 ~ "firebrick3",
      int_value == 4 ~ "red4"
    ))
  } else if (round(ploidy_median) == 3) {
    df <- df %>% dplyr::mutate(colors = case_when(
      int_value == 0 ~ "darkblue",
      int_value == 1 ~ "dodgerblue4",
      int_value == 2 ~ "skyblue2",
      int_value == 3 ~ "azure4",
      int_value == 4 ~ "indianred1",
      int_value == 5 ~ "indianred3",
      int_value == 6 ~ "red4"
    ))
  } else if (round(ploidy_median) == 5) {
    df <- df %>% dplyr::mutate(colors = case_when(
      int_value == 0 ~ "darkblue",
      int_value == 1 ~ "dodgerblue4",
      int_value == 2 ~ "dodgerblue3",
      int_value == 3 ~ "skyblue2",
      int_value == 4 ~ "skyblue1",
      int_value == 5 ~ "azure4",
      int_value == 6 ~ "indianred1",
      int_value == 7 ~ "indianred2",
      int_value == 8 ~ "indianred3",
      int_value == 9 ~ "indianred4",
      int_value == 10 ~ "red4"
    ))
  } else if (round(ploidy_median) == 4) {
    df <- df %>% dplyr::mutate(colors = case_when(
      int_value == 0 ~ "darkblue",
      int_value == 1 ~ "dodgerblue4",
      int_value == 2 ~ "dodgerblue3",
      int_value == 3 ~ "skyblue2",
      int_value == 4 ~ "azure4",
      int_value == 5 ~ "indianred1",
      int_value == 6 ~ "indianred3",
      int_value == 7 ~ "indianred4",
      int_value == 8 ~ "red4"
    ))
  } else if (round(ploidy_median)  > 5) {
    message("Calculated median ploidy is bigger than 5. ratios points will not display colors")
    df <- df %>% dplyr::mutate(colors = "azure4")
  }
  
  # saving the names of the cells so I can iterate later
  cell_names <- df$cell %>% unique()
  
  #chromosomes squares background
  ggchr.back <-
    list(geom_rect(
      data = chrom.rects,
      aes(
        xmin = xstart,
        xmax = xend,
        ymin = -Inf,
        ymax = Inf,
        fill = colors
      ),
      alpha = .2
    ), scale_fill_identity()) 
  
  sec_breaks <- c(0, 0.5e9, 1e9, 1.5e9, 2e9, 2.5e9, 3e9)
  sec_labels <- c(0, 0.5, 1, 1.5, 2, 2.5, 3)
  
  ggaes <- list(
    scale_x_continuous(breaks = xbreaks,
                       labels = gsub("chr", "", chrom.names),
                       position = "top",
                       expand = c(0,0),
                       sec.axis = sec_axis(~., breaks = sec_breaks, labels = sec_labels, name = "Genome Position (Gb)")),
    theme_classic(),
    theme(
      axis.text.x = element_text(angle = 0,
                                 vjust = .5,
                                 size = 15),
      axis.text.y = element_text(size = 15),
      axis.title.y.right = element_text(margin = margin(l = 10)),
      legend.position = "none",
      axis.ticks.x = element_blank(),
      axis.title = element_text(size = 18),
      plot.title = element_text(size = 20)
    )
  )
  
  invisible(parallel::mclapply(seq_along(cell_names), function (x) {
    # just get the mean to calculate the CN, should be 1 though...
    mean_ratios_cell <- df %>% dplyr::filter(cell == cell_names[x]) %>% pull(ratios_rat) %>% mean()
    p <- ggplot(df %>% filter(cell == cell_names[x])) +
      geom_point(aes(abspos, ratios_rat, color = colors),
                 shape = 20,
                 size = 2,
                 alpha = .7) +
      geom_line(aes(abspos, seg_mean), col = "black",
                size = 1.5) +
      ggchr.back +
      ggaes +
      xlab("") +
      ylab("ratios") +
      ggtitle(cell_names[x]) +
      scale_y_continuous(sec.axis = sec_axis(~.*ploidy_median/mean_ratios_cell, breaks = seq(0,max_int_value,1), name  = "Copy Number"), 
                         name = "ratios") +
      scale_color_identity()
    
    cowplot::ggsave(paste(cell_names[x], ".pdf", sep = ""), plot = p, device = "pdf", width = 13, height = 5)
    
  }
  , mc.cores = ncpu)
  )
  
}

ratios_plots(seg_data = dat_seg,
            bin_data = dat)
