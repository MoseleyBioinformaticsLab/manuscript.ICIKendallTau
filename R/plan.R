select_functions = rlang::syms(c("var_select", "pca_select"))
#reference_cor = readRDS(here::here("data/recount_adeno_cor.rds"))
fractions = (c(seq(0.01, 0.05, 0.01), seq(0.1, 0.9, 0.05)))
set.seed(1234)

the_plan <-
  drake_plan(
   #set.seed(1234),
   x = seq(1, 10),
   y = seq(1, 10),
   y2 = seq(10, 1),
   where_na = create_na_indices(20),
   
   positive_kt = compare_positive_kt(x, y, where_na),
   negative_kt = compare_negative_kt(x, y2, where_na),
   positive_local_kt = compare_positive_kt(x, y, where_na, perspective = "local"),
   negative_local_kt = compare_negative_kt(x, y2, where_na, perspective = "local"),
   
   positive_low_kt = compare_positive_kt(x, y, where_na, low_indices = TRUE),
   negative_low_kt = compare_negative_kt(x, y2, where_na, low_indices = TRUE),
   positive_low_local_kt = compare_positive_kt(x, y, where_na, low_indices = TRUE, perspective = "local"),
   negative_low_local_kt = compare_negative_kt(x, y2, where_na, low_indices = TRUE, perspective = "local"),
   
   positive_pearson = compare_positive_pearson(x, y, where_na),
   negative_pearson = compare_negative_pearson(x, y2, where_na),
   
   positive_low_pearson = compare_positive_pearson(x, y, where_na, low_indices = TRUE),
   negative_low_pearson = compare_negative_pearson(x, y2, where_na, low_indices = TRUE),
   
   realistic_sample = create_sample(n = 1000),
   realistic_neg_sample = sort(realistic_sample, decreasing = TRUE),
   realistic_na = create_random_na(),
   realistic_positive_kt = compare_positive_kt(realistic_sample, realistic_sample, realistic_na),
   realistic_negative_kt = compare_negative_kt(realistic_sample, realistic_neg_sample, realistic_na),
   
   realistic_positive_pearson = compare_positive_pearson(realistic_sample, realistic_sample, realistic_na),
   realistic_negative_pearson = compare_negative_pearson(realistic_sample, realistic_neg_sample, realistic_na),
   
   kendall_pearson_comparison = target(
     command = {
       rmarkdown::render(knitr_in("doc/kendall_pearson_comparisons.Rmd"))
       file_out("doc/kendall_pearson_comparisons.pdf")
     }
   ),
   
   transcript_data = readRDS(here::here("data/recount_adeno_counts.rds")),
   transcript_pca = prcomp(t(log1p(transcript_data)), center = TRUE, scale. = FALSE),
   
   select_nonrandom_fraction = target(
      select_function(transcript_pca, transcript_data, fraction = frac_value),
      transform = cross(
         frac_value = !!fractions[c(10, 14)],
         select_function = !!select_functions
      )
   ),
   
   run_nonrandom_fraction = target(
      run_fractional_correlation(select_nonrandom_fraction),
      transform = map(select_nonrandom_fraction)
   ),
   
   select_random_fraction = target(
      select_random(transcript_data, fraction = frac_value),
      transform = map(frac_value = !!fractions)
   ),
   
   run_random_fraction = target(
      run_fractional_correlation(select_random_fraction),
      transform = map(select_random_fraction)
   ),
   
   # compare_nonrandom = target(
   #    compare_fractional_correlation(run_nonrandom_fraction, reference_cor),
   #    transform = map(run_nonrandom_fraction)
   # ),
   # 
   # compare_random = target(
   #    compare_fractional_correlation(run_random_fraction, reference_cor),
   #    transform = map(run_random_fraction)
   # ),
   # 
   improve_runtime = target(
      command = {
         rmarkdown::render(knitr_in("doc/improve_runtime.Rmd"))
         file_out("doc/improve_runtime.pdf")
      }
   )
  

)
