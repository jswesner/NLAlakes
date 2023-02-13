library(tidyverse)
library(brms)
library(tidybayes)

# get the fitted models
mod_d = readRDS(file = "models/mod_d.rds")

# get raw data
d_cell_based_logcentered = readRDS(file = "data/d_cell_based_logcentered.rds")

# 1) set up data grid
newdata = tibble(temp_c = c(-2, 0, 2)) %>%
  expand_grid(pho_c = c(-1, 0 ,1)) %>% 
  expand_grid(zoo_c = c(-1, 0, 1)) %>% 
  expand_grid(biovol_c = seq(min(mod_d$data$biovol_c), max(mod_d$data$biovol_c), length.out = 20))

# 2) sample posteriors
posts_conditional = newdata %>% 
  add_epred_draws(mod_d)

# 3) summarize posteriors
posts_summarized = posts_conditional %>% 
  group_by(temp_c, zoo_c, biovol_c, pho_c) %>% 
  median_qi(.epred)

# 4) plot posteriors
posts_summarized %>% 
  ggplot(aes(x = biovol_c, y = .epred, color = temp_c,
             group = temp_c)) + 
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = temp_c), alpha = 0.5) +
  facet_grid(pho_c~zoo_c)


# 5) plot posteriors again, but with backtransformed values (easier to interpret)
# get means to "uncenter"
means = d_cell_based_logcentered %>% distinct(biovol_mean, phos_mean, zoo_mean, temp_global_mean)

# uncenter and de-log
posts_untransformed = posts_summarized %>% 
  mutate(phos_mean = means$phos_mean,
         zoo_mean = means$zoo_mean,
         biovol_mean = means$biovol_mean,
         temp_mean = means$temp_global_mean) %>% 
  mutate(phos_uncenter = pho_c + phos_mean,
         zoo_uncenter = zoo_c + zoo_mean,
         biovol_uncenter = biovol_c + biovol_mean,
         temp_uncenter = temp_c + temp_mean) %>% 
  mutate(phos = round(10^phos_uncenter, 0),
         zoo = round(10^zoo_uncenter, 0),
         biovol = 10^biovol_uncenter)

# plot again
posts_untransformed %>% 
  ggplot(aes(x = biovol, y = .epred, color = temp_uncenter,
             group = temp_uncenter)) + 
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = temp_uncenter), alpha = 0.5) +
  facet_grid(phos~zoo) + 
  scale_x_log10()



