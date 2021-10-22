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
rmarkdown::render(paste(args[1], "/dashboard/cna_dashboard_zn.Rmd", sep = ""),
                  output_format = "flexdashboard::flex_dashboard",
                  output_file = paste(here(),"/",args[3],"/final_result/",args[2],"_cna_dashboard.html",sep=""), 
                  intermediates_dir =paste(here(),"/",args[3],"/final_result/",sep = ""),
			params = list(
		    outdir = args[3],
		    id=args[2]
                  )
	)

