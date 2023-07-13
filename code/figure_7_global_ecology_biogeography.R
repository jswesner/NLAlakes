# # # MAIN RESULTS OF THE PAPER # # #
# Cell-based CCSR
# Set working directory
# setwd('~/Dropbox/Savina_size_scaling_project/3. analysis/')
# setwd("C:/Users/Thomasm/Dropbox/Work stuff/Projects/Lab group/Savina_size_scaling_project/3. analysis")


library(modelsummary)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(dplyr)
library(ggeffects)
library(interactions)
library(ggthemes)

# Dataset
# Replace Nitrate 0 Values with half of apparent detection limit
# Replace temperature values above 50 with NAs
dat <- read.csv("data/ccsr_cellbased_summaries_2023_06_10.csv") %>%
  select(-c(site_id, uid)) %>%
  filter(phyto_mean_biovol > 0) %>%
  mutate(temp_mean = replace(temp_mean, temp_mean > 50, NA),
         nitrate_mean = log10(replace(nitrate_mean, nitrate_mean == 0, 0.001/2)),
         silicate_mean = log10(silicate_mean),
         zoo_total_biomass = log10(zoo_total_density * zoo_mean_biomass),
         zoo_mean_biomass = log10(zoo_mean_biomass),
         zoo_total_density = log10(zoo_total_density),
         phyto_mean_biovol = log10(phyto_mean_biovol),
         phyto_total_density = log10(phyto_total_density),
         phyto_se_biovol = log10(phyto_se_biovol))
dat2 <- data.frame(scale(dat[,c(1,3:10)], scale = FALSE), phyto_se_biovol = dat$phyto_se_biovol)
dat3 <- dat %>%
  mutate(nitrate_mean = replace(nitrate_mean, nitrate_mean == min(nitrate_mean, na.rm = TRUE), NA))
dat3 <- data.frame(scale(dat3[,c(1,3:10)], scale = FALSE), phyto_se_biovol = dat3$phyto_se_biovol)


## FILTERING OUTLIER LAKE
dat <- dat %>%
  filter(phyto_mean_biovol > 0)

# FITTING PREDICTED INTERACTIONS

mod_final <- lm(phyto_total_density ~ phyto_mean_biovol +
                  temp_mean +
                  nitrate_mean +
                  zoo_total_biomass +
                  phyto_mean_biovol:temp_mean +
                  phyto_mean_biovol:zoo_total_biomass +
                  phyto_mean_biovol:nitrate_mean +
                  temp_mean:zoo_total_biomass +
                  temp_mean:nitrate_mean +
                  #zoo_total_biomass:nitrate_mean +
                  # zoo_total_biomass:nitrate_mean:temp_mean +
                  phyto_mean_biovol:nitrate_mean:temp_mean,
                # zoo_total_biomass:phyto_mean_biovol:temp_mean +
                # zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean,
                data = dat)

summary(mod_final)


# Just for reference: two models with just all 2-way interactions and with all possible interactions
mod_all <- lm(phyto_total_density ~ phyto_mean_biovol +
                temp_mean +
                nitrate_mean +
                zoo_total_biomass +
                phyto_mean_biovol:temp_mean +
                phyto_mean_biovol:zoo_total_biomass +
                phyto_mean_biovol:nitrate_mean +
                temp_mean:zoo_total_biomass +
                temp_mean:nitrate_mean +
                zoo_total_biomass:nitrate_mean +
                zoo_total_biomass:nitrate_mean:temp_mean +
                phyto_mean_biovol:nitrate_mean:temp_mean +
                zoo_total_biomass:phyto_mean_biovol:temp_mean +
                zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean,
              data = dat)


mod_2way <- lm(phyto_total_density ~ phyto_mean_biovol +
                 temp_mean +
                 nitrate_mean +
                 zoo_total_biomass +
                 phyto_mean_biovol:temp_mean +
                 phyto_mean_biovol:zoo_total_biomass +
                 phyto_mean_biovol:nitrate_mean +
                 temp_mean:zoo_total_biomass +
                 temp_mean:nitrate_mean,
               #zoo_total_biomass:nitrate_mean +
               # zoo_total_biomass:nitrate_mean:temp_mean +
               #phyto_mean_biovol:nitrate_mean:temp_mean,
               # zoo_total_biomass:phyto_mean_biovol:temp_mean +
               # zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
               # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean,
               data = dat)

summary(mod_2way)


bbmle::AICtab(mod_2way, mod_final, mod_all)
# mod_final is still best

sjPlot::tab_model(mod_final)

mod_strongpreds <- lm(phyto_total_density ~ phyto_mean_biovol +
                        temp_mean +
                        nitrate_mean +
                        zoo_total_biomass +
                        phyto_mean_biovol:temp_mean +
                        # phyto_mean_biovol:zoo_total_biomass +
                        phyto_mean_biovol:nitrate_mean +
                        # temp_mean:zoo_total_biomass +
                        temp_mean:nitrate_mean +
                        #zoo_total_biomass:nitrate_mean +
                        # zoo_total_biomass:nitrate_mean:temp_mean +
                        phyto_mean_biovol:nitrate_mean:temp_mean,
                      # zoo_total_biomass:phyto_mean_biovol:temp_mean +
                      # zoo_total_biomass:nitrate_mean:phyto_mean_biovol +
                      # phyto_mean_biovol:zoo_total_biomass:nitrate_mean:temp_mean,
                      data = dat)

# library(texPreview)
# library(equatiomatic)
equatiomatic::extract_eq(mod_final)


# mod_final <- lm(phyto_total_density ~ phyto_mean_biovol + temp_mean + nitrate_mean +
#              zoo_total_biomass, data = dat,
#              weights = (1/phyto_se_biovol))


pr <- ggpredict(mod_final, type = "fixed",
                terms = c("phyto_mean_biovol",
                          "temp_mean [8, 20, 32]",
                          # "zoo_total_biomass [1, 2, 3]",
                          "nitrate_mean [-2.4, -1, 0]"))

plot(pr, limits = c(0, 10), show.legend = FALSE)

# # # MAIN FIGURES OF THE PAPER

# # FIGURE 2
# Plotting CCSR Regression Model at cell level
ggplot(data = dat,
       aes(x = phyto_mean_biovol,
           y = phyto_total_density)) + geom_point(color='black',
                                                  alpha = 0.5) + geom_smooth(method = "lm",
                                                                             color='black',
                                                                             se = FALSE) + labs(x = "Log Cell Volume",
                                                                                                y = "Log Density") + theme_few()
# # FIGURE 3
# Plotting Two-Way Interaction Effects of Regression Models

interact_plot(mod_final, pred = phyto_mean_biovol,
              modx = temp_mean,
              mod2= nitrate_mean,
              plot.points = TRUE,
              point.size = 2,
              interval = FALSE,
              modx.values = c(10, 20, 30),
              mod2.values = c(-3, 0),
              colors = c("darkblue", "seagreen", "yellow"),
              x.label = "Log Cell Volume",
              y.label = "Log Density",
              legend.main = "T (°C)") + theme_few() +
  scale_color_viridis(option = "C")

# # FIGURE 4
# Plotting Two-Way Interaction Effects of Regression Models


interact_plot(mod_final, pred = phyto_mean_biovol,
              mod2 = temp_mean,
              modx= nitrate_mean,
              mod2.values = c(10,  30),
              modx.values = c(-3, -1.5, 0),
              colors = c("orange", "purple", "darkblue"),
              plot.points = TRUE,
              point.size = 2,
              interval = FALSE,
              x.label = "Log Cell Volume",
              y.label = "Log Density",
              legend.main = "N (mg/L)") + theme_few() +
  scale_color_viridis(option = "D", direction = -1)

#### Figure 7 ###
newdat <- expand.grid(temp_mean = seq(min(datlow$temp_mean, na.rm = TRUE) - 1 ,
                                      max(datlow$temp_mean, na.rm = TRUE) + 1, 50),
                      phos_total_mean= seq(min(datlow$phos_total_mean, na.rm = TRUE) - 1,
                                           max(datlow$phos_total_mean, na.rm = TRUE) + 1, 50),
                      zoo_total_biomass  = mean(datlow$zoo_total_biomass, na.rm = TRUE),
                      # zoo_total_biomass  = seq(min(dat$zoo_total_biomass, na.rm = TRUE),
                      #                          max(dat$zoo_total_biomass, na.rm = TRUE), 0.1),
                      phyto_mean_biovol  = seq(min(datlow$phyto_mean_biovol, na.rm = TRUE),
                                               max(datlow$phyto_mean_biovol, na.rm = TRUE), 50))

newdat$pred_dens <- predict(mod_final, newdat)
newdat$pred_dens_strongpreds <- predict(mod_strongpreds, newdat)

slopes <- newdat %>%
  group_by(temp_mean, phos_total_mean) %>%
  summarise(slope = coef(lm(pred_dens ~ phyto_mean_biovol))[2])

slopes$slope[slopes$slope > -0.16] <- NA
slopes$slope[slopes$slope < -1.22] <- NA


ggplot(slopes, aes(temp_mean, phos_total_mean)) +
  # geom_raster(aes(fill = slope)) +
  geom_tile(aes(fill = slope)) +
  # scale_fill_viridis_c(option = 'plasma', na.value="grey70") +
  scale_fill_viridis_c(option = 'C', direction = -1, na.value="white") +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(min(datlow$temp_mean, na.rm = TRUE) - 0.5,
                                max(datlow$temp_mean, na.rm = TRUE) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(datlow$phos_total_mean, na.rm = TRUE) - 0.5,
                                max(datlow$phos_total_mean, na.rm = TRUE) + 0.1), expand = c(0, 0)) +
  geom_point(data = dat, col = 'green', alpha = 0.6, size = 2) +
  ylab('Log(Nitrate Concentration)') +
  xlab('Temperature')

## Plotting continues Two-Way Interaction Effects of the slope
dathigh <- dat %>%
  filter(zoo_total_biomass > 2.5)

datlow <- dat %>%
  filter(zoo_total_biomass <= 2.5)

newdat <- expand.grid(temp_mean = seq(min(dathigh$temp_mean, na.rm = TRUE) - 1 ,
                                      max(dathigh$temp_mean, na.rm = TRUE) + 1, 0.1),
                      phos_total_mean= seq(min(dathigh$phos_total_mean, na.rm = TRUE) - 1,
                                           max(dathigh$phos_total_mean, na.rm = TRUE) + 1, 0.02),
                      zoo_total_biomass  = mean(dathigh$zoo_total_biomass, na.rm = TRUE),
                      # zoo_total_biomass  = seq(min(dat$zoo_total_biomass, na.rm = TRUE),
                      #                          max(dat$zoo_total_biomass, na.rm = TRUE), 0.1),
                      phyto_mean_biovol  = seq(min(dathigh$phyto_mean_biovol, na.rm = TRUE),
                                               max(dathigh$phyto_mean_biovol, na.rm = TRUE), 1))

newdat$pred_dens <- predict(mod_final, newdat)
newdat$pred_dens_strongpreds <- predict(mod_strongpreds, newdat)

slopes <- newdat %>%
  group_by(temp_mean, phos_total_mean) %>%
  summarise(slope = coef(lm(pred_dens ~ phyto_mean_biovol))[2])

slopes$slope[slopes$slope > -0.16] <- NA
slopes$slope[slopes$slope < -1.22] <- NA


ggplot(slopes, aes(temp_mean, phos_total_mean)) +
  # geom_raster(aes(fill = slope)) +
  geom_tile(aes(fill = slope)) +
  # scale_fill_viridis_c(option = 'plasma', na.value="grey70") +
  scale_fill_viridis_c(option = 'C', direction = -1, na.value="white") +
  theme_bw(base_size = 18) +
  scale_x_continuous(limits = c(min(dathigh$temp_mean, na.rm = TRUE) - 0.5,
                                max(dathigh$temp_mean, na.rm = TRUE) + 0.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(min(dathigh$phos_total_mean, na.rm = TRUE) - 0.5,
                                max(dathigh$phos_total_mean, na.rm = TRUE) + 0.1), expand = c(0, 0)) +
  geom_point(data = dat, col = 'green', alpha = 0.6, size = 2) +
  ylab('Log(Nitrate Concentration)') +
  xlab('Temperature')

#ggsave('SAR slope ~ temp & nitrate from final model.pdf', height = 8, width = 10)

#
# slopes2 <- newdat %>%
#   group_by(temp_mean, nitrate_mean) %>%
#   summarise(slope = coef(lm(pred_dens_strongpreds ~ phyto_mean_biovol))[2])
#
# slopes2$slope[slopes2$slope > -0.16] <- NA
# slopes2$slope[slopes2$slope < -1.22] <- NA
#
# ggplot(slopes2, aes(temp_mean, nitrate_mean)) +
#   geom_tile(aes(fill = slope)) +
#   scale_fill_viridis_c(option = 'inferno', na.value="grey70") +
#   theme_bw(base_size = 18) +
#   scale_x_continuous(limits = c(min(dat$temp_mean, na.rm = TRUE) - 0.1,
#                                 max(dat$temp_mean, na.rm = TRUE) + 0.1), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(min(dat$nitrate_mean, na.rm = TRUE) - 0.1,
#                                 max(dat$nitrate_mean, na.rm = TRUE) + 0.1), expand = c(0, 0)) +
#   geom_point(data = dat, col = 'green', alpha = 0.6, size = 2) +
#   ylab('log(nitrate concentration)') +
#   xlab('Temperature')
#
# ggsave('SAR slope ~ temp & nitrate from model with only strong predictors.pdf', height = 8, width = 10)



