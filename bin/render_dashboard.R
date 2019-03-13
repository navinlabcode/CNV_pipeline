#!/usr/bin/env Rscript

packages = c("flexdashboard", "here", "rmarkdown", "fs")
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
#  if (!require(x, character.only = TRUE)) {
 #   install.packages(x, dependencies = TRUE)
  #  library(x, character.only = TRUE)
  #}

tryCatch(expr = { library(x,character.only = TRUE)}, 
         error = function(e) { 
         install.packages(x,dependencies=TRUE,type = "source",lib=.libPaths()[1],repos='http://cran.us.r-project.org')}, 
         finally = library(x,character.only = TRUE))

})

args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RSTUDIO_PANDOC= "/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render(paste(args[1], "/dashboard/cna_dashboard.Rmd", sep = ""),
                  output_format = "flexdashboard::flex_dashboard",
                  output_file = paste(dir_ls(here(), 
                                             type = "directory",
                                             regexp = "final_result$", 
                                             recursive = T), 
                                      "/cna_dashboard.html",
                                      sep = ""),
                  params = list(
                    directory = paste(args[1], "/bin", sep = "")
                  ))
