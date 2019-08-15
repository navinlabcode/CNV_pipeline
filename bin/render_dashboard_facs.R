#!/usr/bin/env Rscript

packages = c("flexdashboard", "here", "rmarkdown", "fs")
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {

tryCatch(expr = { library(x,character.only = TRUE)}, 
         error = function(e) { 
         install.packages(x,dependencies=TRUE,type = "source",lib=.libPaths()[1],repos='http://cran.us.r-project.org')}, 
         finally = library(x,character.only = TRUE))

})

args <- commandArgs(trailingOnly = TRUE)

Sys.setenv(RSTUDIO_PANDOC= "/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render(paste(args[1], "/dashboard/cna_dashboard_facs.Rmd", sep = ""),
                  output_format = "flexdashboard::flex_dashboard",
                  output_file = paste(dir_ls(here(), 
                                             type = "directory",
                                             regexp =paste0(args[3],"/final_result$"), 
                                             recurse = T), 
                                      paste0("/",args[2], "_cna_dashboard.html"),
                                      sep = ""),
                  params = list(
                    directory = paste(args[1], "/bin", sep = ""),
		    outdir = args[3],
		    facs = args[4]
                  ))
