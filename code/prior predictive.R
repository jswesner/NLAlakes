library(tidyverse)
library(tidybayes)
library(brms)
library(janitor)
library(ggridges)
library(ggthemes)

# create "plot" folder if one doesn't already exist
if (!dir.exists("plots")){
  dir.create("plots")
} else {
  print("Dir already exists!")
}


# prior simulation --------------------------------------------------------
# load model and data
mod_final_bayessite <- readRDS("models/mod_final_bayessite.rds")
data = mod_final_bayessite$data

biovol_seq = seq(min(data$biovol_c), max(data$biovol_c), length.out = 20)
temp_seq = seq(min(data$temp_c), max(data$temp_c), length.out = 20)

priors = tibble(b = rnorm(100, 0, 1),
               intercept = rnorm(100 ,4, 1),
               sigma = rexp(100, 1),
               biovol = rnorm(100, -0.875, 0.2),
               sims = 1:100)

prior_sims = priors %>%
  expand_grid(biovol_c = biovol_seq) %>%
  expand_grid(temp_c = temp_seq) %>%
  mutate(bio_temp = biovol_c*temp_c) %>%
  mutate(log_phyto = intercept + temp_c*b + biovol_c*biovol + b*temp_c*biovol_c)

prior_sims %>%
  ggplot(aes(x = biovol_c, y = log_phyto)) +
  geom_point()

# load model and data
mod_final_bayessite <- readRDS("models/mod_final_bayessite.rds")
data = mod_final_bayessite$data

mod_final_bayessite = brm(phyto_total_density ~ biovol_c +
                            temp_c +
                            phos_c +
                            zoo_total_c +
                            phos_c +
                            biovol_c:temp_c +
                            biovol_c:zoo_total_c +
                            biovol_c:phos_c +
                            temp_c:zoo_total_c +
                            zoo_total_c:phos_c +
                            biovol_c:phos_c:temp_c+
                            biovol_c:zoo_total_c:temp_c +
                            biovol_c:phos_c:zoo_total_c +
                            temp_c:phos_c:zoo_total_c +
                            temp_c:phos_c:zoo_total_c:biovol_c,
                          family = gaussian(),
                          data = data,
                          prior = c(prior(normal(0, 1), class = "b"),
                                    prior(normal(4, 1), class = "Intercept"),
                                    prior(exponential(1), class = "sigma"),
                                    prior(normal(-0.875, 0.2), coef = "biovol_c")))

saveRDS(mod_final_bayessite, file = "models/mod_final_bayessite.rds")

mod_final_bayessite_priorsens = brm(phyto_total_density ~ biovol_c +
                            temp_c +
                            phos_c +
                            zoo_total_c +
                            phos_c +
                            biovol_c:temp_c +
                            biovol_c:zoo_total_c +
                            biovol_c:phos_c +
                            temp_c:zoo_total_c +
                            zoo_total_c:phos_c +
                            biovol_c:phos_c:temp_c+
                            biovol_c:zoo_total_c:temp_c +
                            biovol_c:phos_c:zoo_total_c +
                            temp_c:phos_c:zoo_total_c +
                            temp_c:phos_c:zoo_total_c:biovol_c,
                          family = gaussian(),
                          data = data,
                          prior = c(prior(normal(0, 2), class = "b"),
                                    prior(normal(0, 5), class = "Intercept"),
                                    prior(exponential(2), class = "sigma"),
                                    prior(normal(0, 2), coef = "biovol_c")))


saveRDS(mod_final_bayessite_priorsens, file = "models/mod_final_bayessite_priorsens.rds")


pp_check(mod_final_bayessite)

mod_final_bayessite = readRDS(file = "models/mod_final_bayessite.rds")
# simulate priors
N = 2000
priors = tibble(b_biovol_c = rnorm(N, -0.875, 0.2),
                b_biovol_c_nitrate_c = rnorm(N, 0, 1),
                b_biovol_c_temp_c = rnorm(N, 0, 1),
                b_biovol_c_temp_c_nitrate_c = rnorm(N, 0, 1),
                b_biovol_c_zoo_total_c = rnorm(N, 0, 1),
                b_nitrate_c = rnorm(N, 0, 1),
                b_phos_c = rnorm(N, 0, 1),
                b_temp_c = rnorm(N, 0, 1),
                b_temp_c_nitrate_c = rnorm(N, 0, 1),
                b_temp_c_zoo_total_c = rnorm(N, 0, 1),
                b_zoo_total_c = rnorm(N, 0, 1),
                b_intercept = rnorm(N, 4, 1),
                sigma = rexp(N, 1)) %>%
  mutate(source = "priors",
         .iter = 1:nrow(.))

# get posteriors
posteriors = mod_final_bayessite %>%
  as_draws_df() %>%
  select(-.iteration, -.chain, -.draw) %>%
  select(!contains(c("lprior", "r_uid"))) %>%
  clean_names() %>%
  select(-lp) %>%
  mutate(source = "posteriors",
         .iter = 1:nrow(.))

# combine priors and posteriors
priors_posteriors = bind_rows(priors, posteriors) %>%
  pivot_longer(cols = c(-source, -.iter))


prior_v_post = priors_posteriors %>%
  group_by(name) %>%
  mutate(median = median(value)) %>%
  ggplot(aes(x = value, y = reorder(name, median), fill = source)) +
  geom_density_ridges(scale = 1) +
  theme_default() +
  scale_fill_grey() +
  labs(fill = "",
       y = "Model Parameter",
       x = "Parameter Value")

ggsave(prior_v_post, file = "plots/prior_v_post.jpg", width = 6, height = 5, dpi = 500)


# prior v posterior for biovolume slope
prior_post_regression = priors_posteriors %>%
  filter(name %in% c("b_biovol_c", "b_intercept")) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  expand_grid(biovol_c = seq(min(data$biovol_c), max(data$biovol_c), length.out = 50)) %>%
  mutate(phyto_abund = b_intercept + b_biovol_c*biovol_c)


prior_post_regression_plot = prior_post_regression %>%
  mutate(source = case_when(source == "priors" ~ "a) prior",
                            TRUE ~ "b) posterior")) %>%
  filter(.iter <= 500) %>%
  ggplot(aes(x = biovol_c, y = phyto_abund)) +
  geom_line(aes(group = .iter), alpha = 0.2) +
  geom_point(data = data %>% mutate(source = "b) posterior"), size = 0.2, aes(y = phyto_total_density)) +
  facet_wrap(~source) +
  theme_default() +
  labs(fill = "",
       y = "log(Phytoplankton abundance)",
       x = "log(Biovolume) centered")

ggsave(prior_post_regression_plot, file = "plots/prior_post_regression_plot.jpg", width = 6, height = 3, dpi = 500)



# prior sensitivity
# get posteriors
posteriors = mod_final_bayessite %>%
  as_draws_df() %>%
  select(-.iteration, -.chain, -.draw) %>%
  select(!contains(c("lprior", "r_uid"))) %>%
  clean_names() %>%
  select(-lp) %>%
  mutate(source = "informative priors",
         .iter = 1:nrow(.))

posteriors_weakpriors = mod_final_bayessite_priorsens %>%
  as_draws_df() %>%
  select(-.iteration, -.chain, -.draw) %>%
  select(!contains(c("lprior", "r_uid"))) %>%
  clean_names() %>%
  select(-lp) %>%
  mutate(source = "weak priors",
         .iter = 1:nrow(.))

library(ggthemes)
prior_sensitivity = bind_rows(posteriors_weakpriors, posteriors) %>%
  pivot_longer(cols = c(-source, -.iter)) %>%
  group_by(name, source) %>%
  median_qi(value) %>%
  ggplot(aes(x = value, y = reorder(name, value), shape = source, color = source)) +
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  position = position_dodge(width = -0.4),
                  size = 0.3) +
  theme_default() +
  # facet_wrap(~source) +
  scale_color_colorblind() +
  # guides(shape = "none") +
  labs(color = "",
       shape = "",
       y = "Model Parameter",
       x = "Parameter Value") +
  theme(legend.position = "top")


ggview::ggview(prior_sensitivity, width = 6, height = 6, units = "in")
ggsave(prior_sensitivity, width = 6, height = 6, units = "in", file = "plots/prior_sensitivity.jpg",
       dpi = 600)

