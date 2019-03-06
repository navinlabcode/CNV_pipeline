packages = c("flexdashboard", "here", "rmarkdown")
# if a package is installed, it will be loaded
# if any are not, the missing package(s) will be installed and loaded
package.check <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

Sys.setenv(RSTUDIO_PANDOC= "/usr/lib/rstudio-server/bin/pandoc")
rmarkdown::render(here("dashboard", "cna_dashboard.Rmd"))
