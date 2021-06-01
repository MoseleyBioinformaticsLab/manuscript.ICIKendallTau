create_complexity_figure = function(single_core_perf){
  pearson = single_core_perf %>%
    dplyr::filter(method %in% "r_pearson")
  
  p_comp = ggplot(pearson, aes(x = n, y = time)) + geom_line() +
    geom_smooth(method = "lm", formula = y ~ (x)) +
    labs(subtitle = "r_pearson, O(n)",
         x = "Number of Features",
         y = "time (s)")
  
  kendall = single_core_perf %>%
    dplyr::filter(method %in% "r_kendall")
  
  k_comp = ggplot(kendall, aes(x = n, y = time)) + geom_line() +
    geom_smooth(method = "lm", formula = y ~ I(x^2)) +
    labs(subtitle = "r_kendall, O(n^2)",
         x = "Number of Features",
         y = "time (s)")
  
  ici = single_core_perf %>%
    dplyr::filter(method %in% "ici_kt")
  
  ici_comp = ggplot(ici, aes(x = n, y = time)) + geom_line() +
    geom_smooth(method = "lm", formula = y ~ x * log(x)) +
    labs(subtitle = "ici_kt, O(nlog(n))",
         x = "Number of Features",
         y = "Time (s)")
  
  p_comp / ici_comp / k_comp
  
}

