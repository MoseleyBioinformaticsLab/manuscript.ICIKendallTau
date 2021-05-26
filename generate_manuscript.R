source("./packages.R")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

rmarkdown::render("./doc/ici_kt_manuscript.Rmd")
