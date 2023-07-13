library(tidyverse)
library(brms)
library(tidybayes)
theme_set(theme_default())

# get the fitted models
mod_d = readRDS(file = "models/mod_d.rds")
mod_nomix = readRDS("models/mod_nomix.rds")
mod_d_colony = readRDS("models/mod_d_colony.rds")
mod_univarate = readRDS("models/mod_univariate.rds")
mod_nomix_univariate = readRDS("models/mod_nomix_univariate.rds")
mod_univariate_colony = readRDS("models/mod_univariate_colony.rds")


# get raw data
d_cell_based_logcentered = readRDS(file = "data/d_cell_based_logcentered.rds") %>%
  mutate(temp_sd = sd(temp_mean, na.rm = T),
         model = "a) Main model")
d_nomix_logcentered = mod_nomix$data %>% mutate(model = "b) No mixotrophs")
d_colony_based_logcentered = readRDS(file = "data/d_colony_based_logcentered.rds") %>%
  mutate(temp_sd = sd(temp_mean, na.rm = T),
         model = "c) Colonial only")

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




# compare colony and no mixotrophs ----------------------------------------
d_draws = as_draws_df(mod_d) %>% mutate(model = "Main model")
nomix_draws = as_draws_df(mod_nomix) %>% mutate(model = "No mixotrophs")
colony_draws = as_draws_df(mod_d_colony) %>% mutate(model = "Colonial taxa only")

all_draws = bind_rows(d_draws, nomix_draws, colony_draws) %>%
  select(-lp__, -.chain, -lprior, -.iteration) %>%
  pivot_longer(cols = c(-.draw, -model))

library(ggridges)
all_draws %>%
  ggplot(aes(x = value, y = name, fill = model)) +
  geom_density_ridges()


biovol = tibble(biovol_c = seq(min(d_cell_based_logcentered$biovol_c),
                      max(d_cell_based_logcentered$biovol_c),
                      length.out = 20))

univ_preds = biovol %>%
  add_epred_draws(mod_univarate) %>%
  group_by(biovol_c) %>%
  median_qi(.epred) %>%
  mutate(log_biovol = biovol_c + unique(d_cell_based_logcentered$biovol_mean),
         model = "a) Main model")

univ_nomix_preds = biovol %>%
  add_epred_draws(mod_nomix_univariate) %>%
  group_by(biovol_c) %>%
  median_qi(.epred) %>%
  mutate(log_biovol = biovol_c + unique(d_cell_based_logcentered$biovol_mean),
         model = "b) No mixotrophs")

univ_colony_preds = biovol %>%
  add_epred_draws(mod_univariate_colony) %>%
  group_by(biovol_c) %>%
  median_qi(.epred) %>%
  mutate(log_biovol = biovol_c + unique(d_cell_based_logcentered$biovol_mean),
         model = "c) Colonial only")

all_preds_uni = bind_rows(univ_preds,
                          univ_nomix_preds,
                          univ_colony_preds)

intercepts = all_preds_uni %>% distinct(model) %>%
  mutate(intercept = c(4.49, 4.49, 3.93))

xvals = tibble(biovol_c= c(-3, -2, -1, 0, 1, 2)) %>%
  mutate(log_biovol = biovol_c + unique(d_cell_based_logcentered$biovol_mean))

raw_data = bind_rows(d_cell_based_logcentered,
                     d_colony_based_logcentered,
                     d_nomix_logcentered)

mod_univariate_colony
mod_nomix_univariate
mod_univarate

slopes = all_preds_uni %>% distinct(model) %>%
  mutate(slope = c(round(summary(mod_univariate)$fixed[[2,1]], 2),
                   round(summary(mod_nomix_univariate)$fixed[[2,1]], 2),
                   round(summary(mod_univariate_colony)$fixed[[2,1]],2)),
         slope_text = paste0("Slope = ", slope))

compare_univariates = all_preds_uni %>%
  ggplot(aes(x = biovol_c, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  facet_wrap(~model) +
  geom_abline(data = intercepts, aes(intercept = intercept,
                                     slope = -0.75),
              linetype = "dotted") +
  geom_point(data = raw_data, aes(y = log_phyto),
             size = 0.2) +
  geom_text(data = slopes, aes(x = -1, y = 0.5, label = slope_text)) +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2),
                     labels = c("-1", "0", "1", "2", "3", "4")) +
  labs(x = "Log10 Cell Biovolume",
       y = "Log10 Density")

ggsave(compare_univariates, file = "plots/compare_univariates.jpg",
       width = 6, height = 2.2)



# make manuscript figures -------------------------------------------------
mean_temp = unique(d_cell_based_logcentered$temp_global_mean)
sd_temp = unique(d_cell_based_logcentered$temp_sd)

d = d_cell_based_logcentered

# prediction grid based on manuscript values in Figures
temp_preds = tibble(temp = c(10, 20, 30)) %>% mutate(temp_c = (temp - mean_temp)/sd_temp)
phos_preds = tibble(pho = c(3, 300)) %>% mutate(log_phos = log10(pho),
                                                pho_c = log_phos - unique(d$phos_mean),
                                                pho_cat = c("Low Res", "High Res"))
zoo_preds = tibble(zoo = c(5, 50)) %>% mutate(log_zoo = log10(zoo),
                                              zoo_c = log_zoo - unique(d$zoo_mean),
                                              zoo_cat = c("Low Pred", "High Pred"))

draws_d = tibble(biovol_c = seq(min(d$biovol_c), max(d$biovol_c), length.out = 20)) %>%
  expand_grid(temp_preds, phos_preds, zoo_preds) %>%
  add_epred_draws(mod_d)

draws_colony = tibble(biovol_c = seq(min(d$biovol_c), max(d$biovol_c), length.out = 20)) %>%
  expand_grid(temp_preds, phos_preds, zoo_preds) %>%
  add_epred_draws(mod_d_colony)

draws_nomix = tibble(biovol_c = seq(min(d$biovol_c), max(d$biovol_c), length.out = 20)) %>%
  expand_grid(temp_preds, phos_preds, zoo_preds) %>%
  add_epred_draws(mod_nomix)

original_interaction = draws_d %>%
  group_by(temp, pho, zoo, biovol_c, pho_cat, zoo_cat) %>%
  median_qi(.epred) %>%
  ggplot(aes(x = biovol_c, y = .epred, color = temp)) +
  geom_line(aes(group = interaction(temp, pho, zoo)),
            size = 1.2) +
  facet_grid(zoo_cat ~ pho_cat) +
  viridis::scale_color_viridis(option = "C")  +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2),
                     labels = c("-1", "0", "1", "2", "3", "4")) +
  labs(x = "Log10 Cell Biovolume",
       y = "Log10 Density",
       color = "Temperature\n(\u00b0 C)",
       title = "a) Main model")

nomix_interaction = draws_nomix %>%
  group_by(temp, pho, zoo, biovol_c, pho_cat, zoo_cat) %>%
  median_qi(.epred) %>%
  ggplot(aes(x = biovol_c, y = .epred, color = temp)) +
  geom_line(aes(group = interaction(temp, pho, zoo)),
            size = 1.2) +
  facet_grid(zoo_cat ~ pho_cat) +
  viridis::scale_color_viridis(option = "C")  +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2),
                     labels = c("-1", "0", "1", "2", "3", "4")) +
  labs(x = "Log10 Cell Biovolume",
       y = "Log10 Density",
       color = "Temperature\n(\u00b0 C)",
       title = "b) No mixotrophs")

colony_interaction = draws_colony %>%
  group_by(temp, pho, zoo, biovol_c, pho_cat, zoo_cat) %>%
  median_qi(.epred) %>%
  ggplot(aes(x = biovol_c, y = .epred, color = temp)) +
  geom_line(aes(group = interaction(temp, pho, zoo)),
            size = 1.2) +
  facet_grid(zoo_cat ~ pho_cat) +
  viridis::scale_color_viridis(option = "C")  +
  scale_x_continuous(breaks = c(-3, -2, -1, 0, 1, 2),
                     labels = c("-1", "0", "1", "2", "3", "4")) +
  labs(x = "Log10 Cell Biovolume",
       y = "Log10 Density",
       color = "Temperature\n(\u00b0 C)",
       title = "c) Colonial only")

method_interaction_compare = original_interaction/nomix_interaction/colony_interaction
ggsave(method_interaction_compare, file = "plots/method_interaction_compare.jpg",
       width = 5, height = 8.5)



# heat map ----------------------------------------------------------------

temp_preds = tibble(temp = seq(min(d$temp_mean, na.rm = T),
                                 max(d$temp_mean, na.rm = T), length.out = 20)) %>% mutate(temp_c = (temp - mean_temp)/sd_temp)
phos_preds = tibble(pho_c = seq(min(d$pho_c, na.rm = T),
                                max(d$pho_c, na.rm = T), length.out = 20)) %>% mutate(log_phos = pho_c - unique(d$phos_mean),
                                                pho = 10^log_phos)
zoo_preds = tibble(zoo = c(5, 50)) %>% mutate(log_zoo = log10(zoo),
                                              zoo_c = log_zoo - unique(d$zoo_mean),
                                              zoo_cat = c("a) Low Pred", "b) High Pred"))

draws_d = tibble(biovol_c = c(0, 1)) %>%
  expand_grid(temp_preds) %>%
  expand_grid(phos_preds) %>%
  expand_grid(zoo_preds) %>%
  add_epred_draws(mod_d) %>% ungroup %>%
  select(-.row, -.chain, -.iteration) %>%
  pivot_wider(names_from = biovol_c, values_from = .epred) %>%
  mutate(slope = `1` - `0`,
         model = "Main model") %>%
  group_by(pho, temp, zoo_cat, model, log_phos, log_zoo, temp_c, pho_c) %>%
  summarize(mean_slopes = mean(slope))

draws_colony = tibble(biovol_c = c(0, 1)) %>%
  expand_grid(temp_preds, phos_preds, zoo_preds) %>%
  add_epred_draws(mod_d_colony) %>% ungroup %>%
  select(-.row, -.chain, -.iteration) %>%
  pivot_wider(names_from = biovol_c, values_from = .epred) %>%
  mutate(slope = `1` - `0`,
         model = "Colonial only") %>%
  group_by(pho, temp, zoo_cat, model, log_phos, log_zoo, temp_c, pho_c) %>%
  summarize(mean_slopes = mean(slope))

draws_nomix = tibble(biovol_c = c(0, 1)) %>%
  expand_grid(temp_preds, phos_preds, zoo_preds) %>%
  add_epred_draws(mod_nomix) %>% ungroup %>%
  select(-.row, -.chain, -.iteration) %>%
  pivot_wider(names_from = biovol_c, values_from = .epred) %>%
  mutate(slope = `1` - `0`,
         model = "No mixotrophs") %>%
  group_by(pho, temp, zoo_cat, model, log_phos, log_zoo, temp_c, pho_c) %>%
  summarize(mean_slopes = mean(slope))

bind_rows(draws_d, draws_colony, draws_nomix) %>%
  # mutate(zoo_cat = as.factor(zoo_cat),
  #        zoo_cat = fct_relevel(zoo_cat, "Low Pred"),
  #        model = as.factor(model),
  #        model = fct_relevel(model, "Main Model", "No mixotrophs")) %>%
  ggplot(aes(x = temp, y = pho_c)) +
  geom_tile(aes(fill = mean_slopes)) +
  scale_fill_viridis_c(option = 'C', direction = -1, na.value="white") +
  facet_grid(model ~ zoo_cat) +
  geom_point(data = dat_grid,
             aes(x = temp_mean, y = pho_c), col = 'green', alpha = 0.6, size = 2)

