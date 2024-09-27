reviewer1_code = function()
{
  withr::local_seed(1713314672)
  ### Generate some data and calculate correlations ###
  correlations <- vector(mode = "list", length = 3)
  names(correlations) <- c("complete", "missing", "imputed")
  correlations <- lapply(correlations, function(x){
    x <- list(pearson = vector("numeric", 1000),
              spearman = vector("numeric", 1000),
              kendall = vector("numeric", 1000),
              ICI_Kt = vector("numeric", 1000))
    return(x)
  })
  #set.seed(1713314672)
  for (m in 1:1000) {
    #m = 1
    x <- 1:1000 + rnorm(1000) * 5000
    y <- 2 * x + 5 + rnorm(1000) * 5000
    df <- data.frame(x = x,
                     y = y)
    for (i in 1:3) {
      #i = 1
      if (i == 1) {
        data <- df
      } else if (i == 2) {
        data <- df
        data$x[data$x < -5000] <- NA
        data$y[data$y < -1000] <- NA
      } else if (i == 3) {
          data <- df
          data$x[data$x < -5000] <- -1e99
          data$y[data$y < -1000] <- -1e99
      }
      for (k in 1:length(correlations[[i]])) {
        if (k %in% c(1, 2)) {
          correlations[[i]][[k]][m] <- cor(na.omit(data)$x, na.omit(data)$y, method =
                                            names(correlations[[i]])[k])
        } else if (k == 3)  {
          # for the cases examined, this should be completely equivalent to method = "kendall", but much faster
          correlations[[i]][[k]][m] <- ici_kt(na.omit(data)$x, na.omit(data)$y)[1]
        } else if (k == 4) {
          correlations[[i]][[k]][m] <- ici_kendalltau(t(data), global_na = NA)$cor['x', 'y']
        }
      }
    }
  }
  return(correlations)
}
