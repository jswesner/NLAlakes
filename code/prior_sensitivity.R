library(tidyverse)
library(tidybayes)
library(brms)
library(janitor)
library(ggridges)
library(ggthemes)

# load fitted models
mod_d = readRDS("models/mod_d.rds")
mod_refit_priorsens = readRDS("models/mod_refit_priorsens.rds")
mod_nomix = readRDS(file = "models/mod_nomix.rds")


posts_d = as_draws_df(mod_d) %>%
  pivot_longer(cols = c(-.iteration, -.chain)) %>%
  mutate(source = "Informative Priors")

posts_sens = as_draws_df(mod_refit_priorsens) %>%
  select(-lp__) %>%
  pivot_longer(cols = c(-.iteration, -.chain, -.draw)) %>%
  mutate(source = "Weak Priors")

posts_nomix = as_draws_df(mod_nomix) %>%
  select(-lp__) %>%
  pivot_longer(cols = c(-.iteration, -.chain, -.draw)) %>%
  mutate(source = "No mixotrophs")

summarized_posts = bind_rows(posts_d, posts_sens, posts_nomix) %>%
  filter(!grepl("lp__", name)) %>%
  filter(!grepl(".draw", name)) %>%
  filter(!grepl("lprior", name)) %>%
  group_by(source, name) %>%
  median_qi(value)


prior_sensitivity = summarized_posts %>%
  ggplot(aes(y = reorder(name, value),
             x = value, xmin = .lower, xmax = .upper, color = source, shape = source)) +
  geom_pointrange(size = 0.4, position = position_dodge(width = -0.4)) +
  theme_default() +
  theme(legend.position = "top") +
  labs(color = "",
       shape = "",
       y = "Parameter Name",
       x = "Parameter Value") +
  scale_color_colorblind()

ggsave(prior_sensitivity, file = "plots/prior_sensitivity.jpg", width = 6, height = 5, units = "in", dpi = 600)
