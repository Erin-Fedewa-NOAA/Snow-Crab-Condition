#Analyze condition metrics in C.opilio using Bayesian multivariate models

#To do: Need to append start date to include Julian day
#Fourth root transform cpue
#Date/location correct condition data with GAM
#How to incorporate lags in this analysis?

# Author: Erin Fedewa
# last updated: 2/27/23

# load ----
library(tidyverse)
library(lubridate)
library(rstan)
library(brms)
library(bayesplot)
library(marginaleffects)
library(emmeans)
library(MARSS)
library(corrplot)
library(factoextra)
library(patchwork)
library(modelr)
library(broom.mixed)
library(pROC)
library(ggthemes)
library(tidybayes)
library(RColorBrewer)
library(knitr)
library(loo)
library(sjPlot)
source("./script/stan_utils.R")

condition_master <- read.csv("./data/total_FA_master.csv")

################################
#data wrangling  
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"),
         maturity != 1,
         Perc_DWT != "NA") %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         region = as.factor(region),
         station = as.factor(gis_station),
         temperature = as.numeric(gear_temperature),
         fourth.root.cpue = as.numeric(cpue^0.25)) -> model.dat 

#Assess collinearity b/w covariates 
model.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            latitude = mean(mid_latitude),
            fourth.root.cpue = mean(fourth.root.cpue)) -> corr.dat

cor(corr.dat[,3:5]) #All < 0.5
corrplot(cor(corr.dat[,3:5]), method = 'number') 

#Distribution of response variable- check family!
model.dat %>%
  ggplot(aes(Perc_DWT)) + 
  geom_histogram()

#############################################
#Prelim model 1: sex, temperature, cpue

mod1_formula <-  bf(Perc_DWT ~ sex + s(temperature, k = 4) + 
                      s(fourth.root.cpue, k=4) + (1 | year/region)) 

mod1 <- brm(mod1_formula,
               data = model.dat,
                family = gaussian,
               cores = 4, chains = 4, iter = 2500,
               save_pars = save_pars(all = TRUE),
               control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod1, file = "./output/mod1.rds")
mod1 <- readRDS("./output/mod1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod1$fit)
neff_lowest(mod1$fit)
rhat_highest(mod1$fit)
summary(mod1)
bayes_R2(mod1)

#Diagnostic Plots
plot(mod1, ask = FALSE)
plot(conditional_smooths(mod1), ask = FALSE)
mcmc_plot(mod1, type = "areas", prob = 0.95)
mcmc_rhat(rhat(mod1)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(mod1, pars = c("b_Intercept", "bs_ssize_1", "bs_sjulian_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(mod1)) #Effective sample size: All ratios > 0.1

##########################

##CPUE condition effect plot
## 95% CI
ce1s_1 <- conditional_effects(mod1 , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mod1 , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mod1 , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = model.dat, aes(x = fourth.root.cpue, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "CPUE", y = "") +
  theme_bw() 

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(mod1 , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mod1 , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mod1 , effect = "temperature", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temperature
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temperature[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temperature[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temperature[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temperature[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = model.dat, aes(x = temperature, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Temperature (C)", y = "") +
  theme_bw() 



