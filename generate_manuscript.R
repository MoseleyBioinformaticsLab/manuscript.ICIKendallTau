library(callr)
supp_materials = r(
  function(){
    source("./packages.R")
    lapply(list.files("./R", full.names = TRUE), source)
    rmarkdown::render(here::here("doc", "supplemental_materials.Rmd"))
  }, show = TRUE
)

supp_tables = r(
  function(){
    source("./packages.R")
    lapply(list.files("./R", full.names = TRUE), source)
    rmarkdown::render(here::here("doc", "supplemental_tables.Rmd"))
  }, show = TRUE
)

manuscript = r(
  function(){
    source("./packages.R")
    lapply(list.files("./R", full.names = TRUE), source)
    rmarkdown::render(here::here("doc", "ici_kt_manuscript.Rmd"))
  }, show = TRUE
)

#beepr::beep(2)
