#Investigate drivers of condition in NBS snow crab using Bayesian multivariate models

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
#data wrangling- NBS dataset  
condition_master %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(lme == "NBS", 
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"), 
         maturity != 1,
         Total_FA_Conc_WWT > 0) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         region = as.factor(sample_region),
         station = as.factor(gis_station),
         temperature = as.numeric(gear_temperature),
         fourth.root.cpue = as.numeric(cpue^0.25),
         fourth.root.invert = as.numeric(total_benthic_cpue^0.25),
         julian = as.numeric(julian)) -> nbs.dat

#Assess collinearity b/w covariates 
nbs.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            fourth.root.cpue = mean(fourth.root.cpue), 
            fourth.root.invert = mean(fourth.root.invert), 
            julian = mean(julian)) -> corr.dat

cor(corr.dat[,3:6]) #All < 0.6
corrplot(cor(corr.dat[,3:6]), method = 'number') 

#Distribution of response variable - choosing a brms family 
nbs.dat %>%
  ggplot(aes(log(Total_FA_Conc_WWT))) + 
  geom_density() #pretty darn left skewed 

#We'll go with truncated Gaussian, just like EBS models 

#############################################
#EBS Models: 
#Model runs not shown here, but group-level effects structure was explored. Due to the high number of stations\
#containing only 1 crab, 1|station and nested 1|region/station models had convergence issues. 
#We'll go with 1|year/region to at least attempt to account for the repeat sampling design

####################################
#MODEL 1 BASE MODEL: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)

nbs_mod1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      (1 | year/region))  

nbs_mod1 <- brm(nbs_mod1_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs_mod1, file = "./output/nbs_mod1.rds")
nbs_mod1 <- readRDS("./output/nbs_mod1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_mod1$fit)
neff_lowest(nbs_mod1$fit)
rhat_highest(nbs_mod1$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs_mod1, ask = FALSE)
plot(conditional_smooths(nbs_mod1), ask = FALSE)
mcmc_plot(nbs_mod1, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs_mod1)) #Effective sample size: All ratios > 0.1
pp_check(nbs_mod1)

summary(nbs_mod1) #credible intervals for spline variance parameters (sds) don't include 0, let's keep smooths
bayes_R2(nbs_mod1) #R2 = 0.2
loo(nbs_mod1) -> a
plot(a)

###########################
#MODEL 2 BASE MODEL + INVERT: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year interaction

nbs_mod2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + (1 | year/region))  

nbs_mod2 <- brm(nbs_mod2_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))


#Save output
saveRDS(nbs_mod2, file = "./output/nbs_mod2.rds")
nbs_mod2 <- readRDS("./output/nbs_mod2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_mod2$fit)
neff_lowest(nbs_mod2$fit)
rhat_highest(nbs_mod2$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs_mod2, ask = FALSE)
plot(conditional_smooths(nbs_mod2), ask = FALSE)
mcmc_plot(nbs_mod2, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs_mod2)) #Effective sample size: All ratios > 0.1
pp_check(nbs_mod2)

summary(nbs_mod2) 
bayes_R2(nbs_mod2) #R2 = 0.29
loo(nbs_mod2) -> b
plot(b)

# model comparison
loo(nbs_mod1, nbs_mod2, moment_match = TRUE) 
#Looks like benthic invert density increases predictive skill 

#######################
#MODEL 3 BASE MODEL + INVERT + CRAB: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year + crab*year interaction

nbs_mod3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + s(fourth.root.cpue, k=4, by = year) +
                      (1 | year/region))  

nbs_mod3 <- brm(nbs_mod3_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))


#Save output
saveRDS(nbs_mod3, file = "./output/nbs_mod3.rds")
nbs_mod3 <- readRDS("./output/nbs_mod3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_mod3$fit)
neff_lowest(nbs_mod3$fit)
rhat_highest(nbs_mod3$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs_mod3, ask = FALSE)
plot(conditional_smooths(nbs_mod3), ask = FALSE)
mcmc_plot(nbs_mod3, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs_mod3)) #Effective sample size: All ratios > 0.1
pp_check(nbs_mod3)

summary(nbs_mod3) 
bayes_R2(nbs_mod3) #R2 = 0.45
loo(nbs_mod3) -> c
plot(c)

# model comparison
loo(nbs_mod1, nbs_mod2, nbs_mod3, moment_match = TRUE) 
#Snow crab CPUE also increases predictive capacity

################################
#MODEL 4 BASE MODEL + INVERT + CRAB + TEMP: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year + crab*year + temp*year interaction 

nbs_mod4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + s(fourth.root.cpue, k=4, by = year) +
                      s(temperature, k=4, by = year) + (1 | year/region))  

nbs_mod4 <- brm(nbs_mod4_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs_mod4, file = "./output/nbs_mod4.rds")
nbs_mod4 <- readRDS("./output/nbs_mod4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_mod4$fit)
neff_lowest(nbs_mod4$fit)
rhat_highest(nbs_mod4$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs_mod4, ask = FALSE)
plot(conditional_smooths(nbs_mod4), ask = FALSE)
mcmc_plot(nbs_mod4, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs_mod4)) #Effective sample size: All ratios > 0.1
pp_check(nbs_mod4)

summary(nbs_mod4) 
bayes_R2(nbs_mod4) #R2 = 0.49
loo(nbs_mod4) -> d
plot(d)

# model comparison
loo(nbs_mod1, nbs_mod2, nbs_mod3, nbs_mod4, moment_match = TRUE)
#So seems that full model has highest predictive capacity 

####################################