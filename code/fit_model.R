library(tidyverse)
library(brms)
library(tidybayes)

# load fitted model
mod_d = readRDS("models/mod_d.rds")
formula = mod_d$formula
prior = mod_d$prior
d = mod_d$data

# refit model
mod_refit = brm(formula = formula,
                data = d,
                family = gaussian(),
                prior = c(prior(normal(4, 1), class = "Intercept"),
                          prior(normal(-0.875, 0.2), coef = "biovol_c"),
                          prior(normal(0, 1), class = "b"),
                          prior(exponential(1), class = "sigma")),
                file = "models/mod_d.rds")

mod_univariate = update(mod_refit, formula = log_phyto ~ biovol_c,
                        newdata = d,
                        prior = c(prior(normal(4, 1), class = "Intercept"),
                                  prior(normal(-0.875, 0.2), coef = "biovol_c"),
                                  prior(exponential(1), class = "sigma")),
                        file = "models/mod_univariate.rds")

# prior sensitivity
mod_refit_priorsens = brm(formula = formula,
                data = d,
                family = gaussian(),
                prior = c(prior(normal(4, 2), class = "Intercept"),
                          prior(normal(-0.875, 0.42), coef = "biovol_c"),
                          prior(normal(0, 2), class = "b"),
                          prior(exponential(1), class = "sigma")),
                file = "models/mod_refit_priorsens.rds")

saveRDS(mod_refit_priorsens, file = "models/mod_refit_priorsens.rds")



# re-run without mixotrophs -----------------------------------------------
# 1) get temps from older dataset (can we just use the new data set for this?)
d_temp_c = read_csv("ccsr_colonybased_summaries_2022_10_27.csv") %>%
  select(uid, site_id, temp_mean) %>%
  mutate(temp_c = scale(temp_mean))

# load dataset without mixotrophs
d_no_mix = read_csv("data/ccsr_cellbased_summaries_2023_06_10.csv") %>%
  select(-temp_mean) %>%
  left_join(d_temp_c) %>%
  mutate(zoo_total_biomass = zoo_mean_biomass*zoo_total_density)

# transform and center
d_no_mix_logcentered = d_no_mix %>%
  mutate(log_phos = log10(phos_total_mean),
         log_zoo = log10(zoo_total_biomass),
         log_biovol = log10(phyto_mean_biovol),
         log_phyto = log10(phyto_total_density)) %>%
  mutate(phos_mean = mean(log_phos, na.rm = T),
         zoo_mean = mean(log_zoo, na.rm = T),
         biovol_mean = mean(log_biovol, na.rm = T),
         temp_global_mean = mean(temp_mean, na.rm = T)) %>%
  mutate(pho_c = log_phos - phos_mean,
         zoo_c = log_zoo - zoo_mean,
         biovol_c = log_biovol - biovol_mean)

saveRDS(d_no_mix_logcentered, file = "data/d_no_mix_logcentered.rds")

# re-fit model
mod_nomix = readRDS("models/mod_nomix.rds")

mod_nomix = update(mod_d, newdata = d_no_mix_logcentered,
                   file = "models/mod_nomix.rds")

mod_nomix_univariate = brm(formula = log_phyto ~ biovol_c,
                           data = d_no_mix_logcentered,
                           family = gaussian(),
                           prior = c(prior(normal(4, 2), class = "Intercept"),
                                     prior(normal(-0.875, 0.42), coef = "biovol_c"),
                                     prior(exponential(1), class = "sigma")),
                           file = "models/mod_nomix_univariate.rds")


# colony based ----------------------------------------------
mod_d = readRDS("models/mod_d.rds") # load original model to refit with colony data later

# wrangle colony data
d_temp_c = read_csv("ccsr_colonybased_summaries_2022_10_27.csv") %>%
  select(uid, site_id, temp_mean) %>%
  mutate(temp_c = scale(temp_mean))

d_colony_based = read_csv("ccsr_colonybased_summaries_2022_10_27.csv") %>%
  select(-temp_mean) %>%
  left_join(d_temp_c) %>%
  mutate(zoo_total_biomass = zoo_mean_biomass*zoo_total_density)

d_colony_based_logcentered = d_colony_based %>%
  mutate(log_phos = log10(phos_total_mean),
         log_zoo = log10(zoo_total_biomass),
         log_biovol = log10(phyto_mean_biovol),
         log_phyto = log10(phyto_total_density)) %>%
  mutate(phos_mean = mean(log_phos, na.rm = T),
         zoo_mean = mean(log_zoo, na.rm = T),
         biovol_mean = mean(log_biovol, na.rm = T),
         temp_global_mean = mean(temp_mean, na.rm = T)) %>%
  mutate(pho_c = log_phos - phos_mean,
         zoo_c = log_zoo - zoo_mean,
         biovol_c = log_biovol - biovol_mean)

saveRDS(d_colony_based_logcentered, file = "data/d_colony_based_logcentered.rds")


mod_d_colony = update(mod_d, newdata = d_colony_based_logcentered,
                      file = "models/mod_d_colony.rds")

mod_univariate_colony = update(mod_nomix_univariate, newdata = d_colony_based_logcentered,
                               file = "models/mod_univariate_colony.rds")

