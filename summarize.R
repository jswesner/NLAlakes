library(tidyverse)
library(brms)
library(tidybayes)


# univariate models -------------------------------------------------------
mod_univarate = readRDS("models/mod_univariate.rds")
mod_nomix_univariate = readRDS("models/mod_nomix_univariate.rds")
mod_univariate_colony = readRDS("models/mod_univariate_colony.rds")

draws_univarate = as_draws_df(mod_univariate) %>% mutate(model = "Main model")
draws_nomix_univariate = as_draws_df(mod_nomix_univariate) %>% mutate(model = "No mixotrophs")
draws_univarate_colonial = as_draws_df(mod_univariate_colony) %>% mutate(model = "Colonial taxa only")

all_draws_uni = bind_rows(draws_univarate,
                          draws_nomix_univariate,
                          draws_univarate_colonial) %>%
  select(-lp__, -.chain, -lprior, -.iteration) %>%
  pivot_longer(cols = c(-.draw, -model))

all_draws_uni %>%
  group_by(model, name) %>%
  median_qi(value) %>%
  select(-.width, -.point, -.interval) %>%
  mutate(median_qi = paste0(round(value, 1), " (", round(.lower,1), " to ", round(.upper,1), ")")) %>%
  select(-value, -.lower, -.upper) %>%
  pivot_wider(names_from = model, values_from = median_qi) %>%
  write_csv("tables/compare_methods_univariate.csv")


# interaction models ------------------------------------------------------

# load fitted models
mod_d = readRDS(file = "models/mod_d.rds")
mod_nomix = readRDS("models/mod_nomix.rds")
mod_d_colony = readRDS("models/mod_d_colony.rds")

# get draws
d_draws = as_draws_df(mod_d) %>% mutate(model = "Main model")
nomix_draws = as_draws_df(mod_nomix) %>% mutate(model = "No mixotrophs")
colony_draws = as_draws_df(mod_d_colony) %>% mutate(model = "Colonial taxa only")
all_draws = bind_rows(d_draws, nomix_draws, colony_draws) %>%
  select(-lp__, -.chain, -lprior, -.iteration) %>%
  pivot_longer(cols = c(-.draw, -model))

# summarize draws
all_draws %>%
  group_by(model, name) %>%
  median_qi(value) %>%
  select(-.width, -.point, -.interval) %>%
  mutate(median_qi = paste0(round(value, 1), " (", round(.lower,1), " to ", round(.upper,1), ")")) %>%
  select(-value, -.lower, -.upper) %>%
  pivot_wider(names_from = model, values_from = median_qi) %>%
  write_csv("tables/compare_methods.csv")
