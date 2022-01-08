# for the manuscript
source("./packages.R")
lapply(list.files("./R", full.names = TRUE), source)
# we process the supplement firs so we can refer to the figures and tables
# in there in the main manuscript
rmarkdown::render("./doc/supplemental_materials.Rmd")
#beepr::beep(2)
rmarkdown::render("./doc/ici_kt_manuscript.Rmd")

beepr::beep(2)
