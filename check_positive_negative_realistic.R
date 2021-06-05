source("packages.R")
lapply(list.files("./R", full.names = TRUE), source)

loadd(realistic_sample_1)
loadd(realistic_sample_2)
loadd(realistic_neg_sample)

loadd(realistic_na)

n_na = purrr::map_int(realistic_na, ~ length(.x))

sample_df = data.frame(
  s1 = realistic_sample_1,
  s2 = realistic_sample_2,
  neg = realistic_neg_sample
)

ggplot(sample_df, aes(x = s1, y = s2)) + 
  geom_point()
ggplot(sample_df, aes(x = s1, y = neg)) +
  geom_point()

test_realistic_na = realistic_na[(n_na >= 500) & (n_na <= 750)]
test_na = test_realistic_na[[1]]
s1_na = test_na[test_na <= 1000]
s2_na = test_na[test_na > 1000] - 1000
miss_df = sample_df
miss_df$s1[s1_na] = NA
miss_df$s2[s2_na] = NA
miss_df$neg[s2_na] = NA
ggplot(miss_df, aes(x = s1, y = s2)) +
  geom_point()
ggplot(miss_df, aes(x = s1, y = neg)) + 
  geom_point()

ici_kt(miss_df$s1, miss_df$s2, "global")
ici_kt(miss_df$s1, miss_df$neg, "global")

miss_df2 = miss_df
miss_df2[is.na(miss_df$s1), "s1"] = 0.43
miss_df2[is.na(miss_df$neg), "neg"] = -13.88

ggplot(miss_df2, aes(x = s1, y = neg)) +
  geom_point()
ggsave("check_negative_kendall.png", device = "png", width = 8, height = 5, units = "in")

# ok, that's the problem, how can we solve it?
# one solution is to check the direction of correlation,
# and change one of the inputs to use the max or min, instead
# of both using the min
ici_kt_dir = function(x, y, perspective = "global"){
  na_x_y = is.na(x) | is.na(y)
  
  x_no_na = x[!na_x_y]
  y_no_na = y[!na_x_y]
  
  init_cor = ici_kt(x_no_na, y_no_na, perspective = perspective)[1]
  if (init_cor > 0) {
    replace_x = min(x, na.rm = TRUE) - 0.1
    replace_y = min(y, na.rm = TRUE) - 0.1
  } else {
    replace_x = min(x, na.rm = TRUE) - 0.1
    replace_y = max(y, na.rm = TRUE) + 0.1
  }
  
  x_new = x
  x_new[is.na(x)] = replace_x
  
  y_new = y
  y_new[is.na(y)] = replace_y
  
  tmp = ici_kt(x_new, y_new, perspective = perspective)
  list(cor = tmp, data = data.frame(x = x_new, y = y_new))
}

test_1 = ici_kt_dir(miss_df$s1, miss_df$s2)
test_2 = ici_kt_dir(miss_df$s1, miss_df$neg)
