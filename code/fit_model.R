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
                          prior(exponential(1), class = "sigma")))


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
