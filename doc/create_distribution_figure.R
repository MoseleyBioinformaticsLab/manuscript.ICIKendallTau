library(ggplot2)
library(dplyr)
theme_set(theme_void())

set.seed(1234)

d_vals = seq(-5, 5, 0.01)
r_norm = data.frame(d = dnorm(d_vals),
                    x = d_vals) %>%
  mutate(c = cumsum(d))

ggplot(r_norm, aes(x = x, y = d)) + geom_line()
ggplot(r_norm, aes(x = x, y = c)) + geom_line()

r_norm = dplyr::mutate(r_norm, c = cumsum(d))

p = ggplot(r_norm, aes(x = x, y = c)) + geom_line()
p
ggsave(here::here("doc/instrument_response.svg"),
       p, device = "svg")

dl_vals = seq(0, 10, 0.01)
l_norm = data.frame(x = dl_vals,
                    y = dlnorm(dl_vals, meanlog = 1.5, sdlog = 1))

dl = ggplot(l_norm, aes(x = x, y = y)) + geom_line()
dl
ggsave(here::here("doc/sample_values.svg"),
       dl, device = "svg")
