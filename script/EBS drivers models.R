#Investigate drivers of condition in EBS snow crab using Bayesian multivariate models

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
  #2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
  #used as response variable in all further analyses. 

# Author: Erin Fedewa
# last updated: 3/12/24

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
#look at distributions of cpue response variables 
condition_master %>%
  ggplot() +
  geom_histogram(aes(cpue))

condition_master %>%
  ggplot() +
  geom_histogram(aes(total_benthic_cpue))

#data wrangling- EBS dataset  
condition_master %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(lme == "EBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"), 
         maturity != 1) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         region = as.factor(sample_region),
         station = as.factor(gis_station),
         temperature = as.numeric(gear_temperature),
         fourth.root.cpue = as.numeric(cpue^0.25),
         fourth.root.invert = as.numeric(total_benthic_cpue^0.25)) -> ebs.dat 

#Assess collinearity b/w covariates 
ebs.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            latitude = mean(mid_latitude),
            fourth.root.cpue = mean(fourth.root.cpue), 
            fourth.root.invert = mean(fourth.root.invert)) -> corr.dat

cor(corr.dat[,3:6]) #All < 0.6
corrplot(cor(corr.dat[,3:6]), method = 'number') 

#Distribution of response variable - choosing a brms family 
ebs.dat %>%
  ggplot(aes(Total_FA_Conc_WWT)) + 
  geom_density() #pretty darn left skewed 

#Test glm model to look at distribution of residuals 
test.1 <- glm(Total_FA_Conc_WWT ~ temperature, data=ebs.dat, family =gaussian(link=log))
plot(test.1) #overdispersion in qqplot, not gaussian! 
plot(density(resid(test.1, type='deviance'))) #very long tail, much heavier than Gaussian

# fitting a test brms model with a Gaussian likelihood - truncating response at 0, else models predict values < 0
model_normal <- brm(Total_FA_Conc_WWT | trunc(lb = 0) ~ 1, family = gaussian(link="identity"), data = ebs.dat)
summary(model_normal)

# fitting a test brms model with a skew normal likelihood
model_skew <- brm(Total_FA_Conc_WWT ~ 1, family = skew_normal(link="log"), data = ebs.dat)
summary(model_skew)

# fitting a test brms model with a gamma likelihood
model_gamma <- brm(Total_FA_Conc_WWT ~ 1, family = "gamma", data = ebs.dat)
summary(model_gamma)

# fitting a test brms model with a lognormal likelihood
model_log <- brm(Total_FA_Conc_WWT ~ 1, family = lognormal(), data = ebs.dat)
summary(model_log)

# posterior predictive checking
pp_check(model_normal, ndraws = 1e2) + pp_check(model_skew, ndraws = 1e2) +
 pp_check(model_gamma, ndraws = 1e2) + pp_check(model_log, ndraws = 1e2)

# posterior predictive checking - boxplots
pp_check(model_normal, type = "boxplot", ndraws = 20) + pp_check(model_skew, type = "boxplot", ndraws = 20) +
  pp_check(model_gamma, type = "boxplot", ndraws = 20) + pp_check(model_log, type = "boxplot", ndraws = 20)

#let's look at the distribution of minimum values for posterior distributions vrs data
pp_check(model_normal, type = "stat", stat = "min") + pp_check(model_skew, type = "stat", stat = "min") +
  pp_check(model_gamma, type = "stat", stat = "min") + pp_check(model_log, type = "stat", stat = "min")
#skew normal is predicting negative values 

#now maximum values
pp_check(model_normal, type = "stat", stat = "max") + pp_check(model_skew, type = "stat", stat = "max") +
    pp_check(model_gamma, type = "stat", stat = "max") + pp_check(model_log, type = "stat", stat = "max")
#lognormal is way overshooting max 

#and means
pp_check(model_normal, type = "stat", stat = "mean") + pp_check(model_skew, type = "stat", stat = "mean") +
  pp_check(model_gamma, type = "stat", stat = "mean") + pp_check(model_log, type = "stat", stat = "mean")
  
model_normal <- add_criterion(model_normal, "waic")
model_skew <- add_criterion(model_skew, "waic")
model_gamma <- add_criterion(model_gamma, "waic")
model_log <- add_criterion(model_log, "waic")
  loo_compare(model_normal, model_skew, model_gamma, model_log, criterion = "waic")
#predictive accuracy highest for gaussian truncated model
  
#to do- how to add bounds to skew-normal??? 

#############################################
#EBS Models: 1) sqrt with gaussian/"identity link

#Most complex should be - or test density/temp/invert all seperate? 
mod1_formula <-  bf(Perc_DWT ~ sex + cw + s(temperature, k = 4) + s(julian, k = 4)
                      s(fourth.root.cpue, k=4) + s(fourth.root.invert, k=4) + (1 | year/region)) 





#Model 1: %dWT ~ sex, cw, temperature, cpue, benthic invert,  

mod1_formula <-  bf(Perc_DWT ~ sex + cw + s(temperature, k = 4) + 
                      s(fourth.root.cpue, k=4) + s(fourth.root.invert, k=4) + (1 | year/region)) 

mod1 <- brm(mod1_formula,
               data = ebs.dat,
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

pp_check(mod1)

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
  geom_point(data = ebs.dat, aes(x = fourth.root.cpue, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
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
  geom_point(data = ebs.dat, aes(x = temperature, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Temperature (C)", y = "") +
  theme_bw() 

##Benthic Inverts 
## 95% CI
ce1s_1 <- conditional_effects(mod1 , effect = "fourth.root.invert", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mod1 , effect = "fourth.root.invert", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mod1 , effect = "fourth.root.invert", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.invert
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.invert[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.invert[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.invert[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.invert[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = ebs.dat, aes(x = fourth.root.invert, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Benthic Invertebrate Density", y = "") +
  theme_bw() 

#Size effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
## 95% CI
ce1s_1 <- conditional_effects(mod1 , effect = "cw", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mod1 , effect = "cw", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mod1 , effect = "cw", re_formula = NA,
                              probs = c(0.1, 0.9))

dat_ce <- ce1s_1$cw
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$cw[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$cw[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$cw[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$cw[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = ebs.dat, aes(x = cw, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Carapace Width (mm)", y = "") +
  theme_bw() 

#Sex effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
ce1s_1 <- conditional_effects(mod1, effect = "sex", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$sex %>%
  dplyr::select(sex, estimate__, lower__, upper__) %>%
  mutate(sex = case_when(sex == 1 ~ "Male",
                         sex == 2 ~ "Female")) %>%
  ggplot(aes(factor(sex, levels = c("Male", "Female")), estimate__)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), width=0.3, size=0.5) +
  labs(y="% Hepatopancreas Dry Weight", x="") +
  theme_bw() 

########################
#Simpler model- pulling temp b/c that is essentially a year effect 

mod2_formula <-  bf(Perc_DWT ~ sex +  
                      s(fourth.root.cpue, k=4)  + year + (1 | region)) 

mod2 <- brm(mod2_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod2, file = "./output/mod2.rds")
mod2 <- readRDS("./output/mod2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod2$fit)
neff_lowest(mod1$fit)
rhat_highest(mod1$fit)
summary(mod2)
bayes_R2(mod1)

conditional_effects(mod2)
#condition better at high density - better habitat? 

##CPUE condition effect plot
## 95% CI
ce1s_1 <- conditional_effects(mod2 , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(mod2 , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(mod2 , effect = "fourth.root.cpue", re_formula = NA,
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
  geom_point(data = ebs.dat, aes(x = fourth.root.cpue, y = Perc_DWT), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "CPUE", y = "") +
  theme_bw() 

#################################

mod3_formula <-  bf(Perc_DWT ~ sex + year + 
                        year/fourth.root.cpue + (1 | region)) 

mod3 <- brm(mod3_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod3, file = "./output/mod3.rds")
mod3 <- readRDS("./output/mod3.rds")

summary(mod3)
bayes_R2(mod3)

#Diagnostic Plots
plot(mod3, ask = FALSE)
plot(conditional_smooths(mod3))
conditional_effects(mod3)

#pre-collapse, low condition regardless of density, 
#different effect of density pre-collapse than post 

##############################################

mod4_formula <-  bf(Perc_DWT ~ sex +  
                      s(temperature, k=4) + (1 | region)) 

mod4 <- brm(mod4_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod4, file = "./output/mod4.rds")
mod4 <- readRDS("./output/mod4.rds")
summary(mod4)
conditional_effects(mod4)

#temperature is an annual effect - 
#always keep sex in as nuisance variable, region for spatial correlation

#############################################

mod5_formula <-  bf(Perc_DWT ~ sex + s(temperature, k=4) +
                      year/fourth.root.cpue + (1 | region)) 
#try s(fourth.root.cpue, k=4, by=year)

mod5 <- brm(mod5_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod5, file = "./output/mod5.rds")
mod5 <- readRDS("./output/mod5.rds")
summary(mod5)
conditional_effects(mod5)

#Not estimating temperature effect in first model- low condition with different
  #response to density. now to just look at model with temperature b/c
#Cant distinguish between density and temp in same model 


