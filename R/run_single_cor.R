##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

##' @return
##' @author rmflight
##' @export
run_single_cor <- function() {

  # running a bunch of timings on a single cor
  future::plan(multicore(workers = 1))
  
  n_to_run = seq(100, 5000, 100)
  
  r_kendall = function(x, y){
    s_start = Sys.time()
    cor_val = cor(x, y, method = "kendall")
    s_end = Sys.time()
    data.frame(cor = cor_val,
               method = "r_kendall",
               n = length(x),
               time = as.numeric(difftime(s_end, s_start, units = "secs"))
              )
  }
  r_pearson = function(x, y){
    s_start = Sys.time()
    cor_val = cor(x, y, method = "pearson")
    s_end = Sys.time()
    data.frame(cor = cor_val,
               method = "r_pearson",
               n = length(x),
               time = as.numeric(difftime(s_end, s_start, units = "secs"))
    )
  }
  
  i_kendall = function(x, y){
    s_start = Sys.time()
    cor_val = ici_kt(x, y, "global")
    s_end = Sys.time()
    data.frame(cor = cor_val[1],
               method = "ici_kt",
               n = length(x),
               time = as.numeric(difftime(s_end, s_start, units = "secs"))
    )
  }
  
  run_all = purrr::map_df(
    n_to_run, function(n){
      x = rnorm(n)
      y = rnorm(n)
      ici = i_kendall(x, y)
      pear = r_pearson(x, y)
      kend = r_kendall(x, y)
      bind_rows(ici, pear, kend)
    }
  )
  
  future::plan(multicore)

  run_all
}
