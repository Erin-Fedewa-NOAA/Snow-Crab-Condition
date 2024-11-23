#Investigate the effects of temperature and snow crab density on 
  #energetic condition in NBS

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
#2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
#used as response variable in all further analyses. 

# Author: EJF

#load 
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

#--------------------------------------------------------------------------------
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
         cpue = as.numeric(cpue),
         invert = as.numeric(total_benthic_cpue),
         julian = as.numeric(julian)) -> nbs.dat

#Assess collinearity b/w covariates 
nbs.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            cpue = mean(cpue), 
            invert = mean(invert), 
            julian = mean(julian)) -> corr.dat

cor(corr.dat[,3:6]) 
corrplot(cor(corr.dat[,3:6]), method = 'number') 

#Distribution of response variable - choosing a brms family 
nbs.dat %>%
  ggplot(aes(Total_FA_Conc_WWT)) + 
  geom_density() #long tail, but more normalish than ebs dataset  

#We'll go with truncated Gaussian, just like EBS models 

#-----------------------------------------------------------------------
#NBS models: Using same approach as EBS drivers models.R script, but we're not 
  #testing benthic invert models 

#MODEL 1 base model: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)

nbs1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      (1 | region)) 

## Show default priors
get_prior(nbs1_formula, nbs.dat)

nbs1 <- brm(nbs1_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs1, file = "./output/nbs1.rds")
nbs1 <- readRDS("./output/nbs1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs1$fit)
neff_lowest(nbs1$fit)
rhat_highest(nbs1$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs1, ask = FALSE)
plot(conditional_effects(nbs1), ask = FALSE)
mcmc_plot(nbs1, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs1)) #Effective sample size: All ratios > 0.1
pp_check(nbs1) #not bad!

summary(nbs1) #credible intervals for spline variance parameters (sds) don't include 0, let's keep smooths
bayes_R2(nbs1) #R2 = 0.05
loo(nbs1)

#----------------------------------------------------------------------------------
#MODEL 2 BASE + CRAB CPUE: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + crab cpue main effect 

nbs2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(fourth.root.cpue, k = 3) + (1 | region))  

nbs2 <- brm(nbs2_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))


#Save output
saveRDS(nbs2, file = "./output/nbs2.rds")
nbs2 <- readRDS("./output/nbs2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs2$fit)
neff_lowest(nbs2$fit)
rhat_highest(nbs2$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs2, ask = FALSE)
plot(conditional_effects(nbs2), ask = FALSE)
mcmc_plot(nbs2, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs2)) #Effective sample size: All ratios > 0.1
pp_check(nbs2)

summary(nbs2) 
bayes_R2(nbs2) #R2 = 0.06
loo(nbs2) 

# model comparison
loo(nbs1, nbs2, moment_match = TRUE) 
#Really no difference in these two models. We'll keep cpue in as a way to test hypotheses
  #in full model 

#--------------------------------------------------------------------------------
#MODEL 3 BASE + TEMP: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature fixed effect 

nbs3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(temperature, k=3) + (1 | region))  

nbs3 <- brm(nbs3_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs3, file = "./output/nbs3.rds")
nbs3 <- readRDS("./output/nbs3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs3$fit)
neff_lowest(nbs3$fit)
rhat_highest(nbs3$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs3, ask = FALSE)
plot(conditional_effects(nbs3), ask = FALSE)
mcmc_plot(nbs3, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs3)) #Effective sample size: All ratios > 0.1
pp_check(nbs3)

summary(nbs3) 
bayes_R2(nbs3) #R2 = 0.10
loo(nbs3) -> d
plot(d)

# model comparison
loo(nbs1, nbs2, nbs3, moment_match = TRUE)
#Temperature model has highest predictive capacity 

#------------------------------------------------------------------------------
#MODEL 4 BASE + TEMP + DENSITY: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature and density fixed effect 

nbs4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(temperature, k=3) + s(fourth.root.cpue, k=3) +
                      (1 | region))  

nbs4 <- brm(nbs4_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs4, file = "./output/nbs4.rds")
nbs4 <- readRDS("./output/nbs4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs4$fit)
neff_lowest(nbs4$fit)
rhat_highest(nbs4$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs4, ask = FALSE)
plot(conditional_effects(nbs4), ask = FALSE)
mcmc_plot(nbs4, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs4)) #Effective sample size: All ratios > 0.1
pp_check(nbs4)

summary(nbs4) 
bayes_R2(nbs4) #R2 = 0.11
loo(nbs4) 

# model comparison
loo(nbs1, nbs2, nbs3, nbs4, moment_match = TRUE)
#So seems that full additive nbs4 has highest predictive capacity 
#and is a substantial improvement over nbs1 and nbs2 but addition of cpue 
#doesnt improve much 
  
#------------------------------------------------------------------------------
#MODEL 5 BASE + TEMP x DENSITY INTERACTION: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature x density linear interaction

nbs5_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      fourth.root.cpue*temperature + (1 | region))   

nbs5 <- brm(nbs5_formula,
            data = nbs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(nbs5, file = "./output/nbs5.rds")
nbs5 <- readRDS("./output/nbs5.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs5$fit)
neff_lowest(nbs5$fit)
rhat_highest(nbs5$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(nbs5, ask = FALSE)
plot(conditional_effects(nbs5), ask = FALSE)
mcmc_plot(nbs5, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(nbs5)) #Effective sample size: All ratios > 0.1
pp_check(nbs5)

summary(nbs5) 
bayes_R2(nbs5) #R2 = 0.06
loo(nbs5) 

# model comparison
loo(nbs1, nbs2, nbs3, nbs4, nbs5, moment_match = TRUE)
#Looks like interaction did not improve predictive capacity from main effects model

#--------------------------------------------------------------------------------
#Full model Comparison (base model + base/cpue + base/cpue/temp + base/cpue*temp)

#LOO-CV
nbs1 <- add_criterion(nbs1, "loo")
nbs3 <- add_criterion(nbs3, "loo")
nbs4 <- add_criterion(nbs4, "loo")
nbs5 <- add_criterion(nbs5, "loo")
loo_compare(nbs1, nbs3, nbs4, nbs5, criterion = "loo") %>% print(simplify = F)
model.comp <- loo(nbs1, nbs3, nbs4, nbs5, moment_match = TRUE)

#and loo weights
model_weights(nbs1, nbs3, nbs4, nbs5, weights = "loo") %>% round(digits = 2)

#Table of Rsq Values 
rbind(bayes_R2(nbs1), 
      bayes_R2(nbs3),
      bayes_R2(nbs4),
      bayes_R2(nbs5)) %>%
  as_tibble() %>%
  mutate(model = c("nbs1", "nbs3", "nbs4", "nbs5"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#weights 
loo1 <- loo(nbs1)
loo3 <- loo(nbs3)
loo4 <- loo(nbs4)
loo5 <- loo(nbs5)

loo_list <- list(loo1, loo3, loo4, loo5)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
#Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)

#Save model output 
tab_model(nbs1, nbs2, nbs3, nbs4, nbs5)

forms <- data.frame(formula=c(as.character(nbs1_formula)[1],
                              as.character(nbs2_formula)[1],
                              as.character(nbs3_formula)[1],
                              as.character(nbs4_formula)[1],
                              as.character(nbs5_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/nbs_pop_model_comp.csv")

#--------------------------------------------------------------------------
#FINAL MODEL:  Run nbs4 model with 10,000 iterations and set seed for reproducibility 
nbs_pop_final <- brm(nbs4_formula,
                     data = nbs.dat,
                     family = gaussian,
                     cores = 4, chains = 4, iter = 10000, warmup = 1000,
                     save_pars = save_pars(all = TRUE), seed = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(nbs_pop_final, file = "./output/nbs_pop_final.rds")
nbs_pop_final <- readRDS("./output/nbs_pop_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_pop_final$fit)
neff_lowest(nbs_pop_final$fit)
rhat_highest(nbs_pop_final$fit)
summary(nbs_pop_final)
bayes_R2(nbs_pop_final) #r2 = .11
loo1 <- loo(nbs_pop_final, save_psis = TRUE)
plot(loo1)

#Diagnostic Plots
plot(nbs_pop_final, ask = FALSE)
plot(conditional_effects(nbs_pop_final), ask = FALSE)
mcmc_plot(nbs_pop_final, prob = 0.95)
mcmc_neff(neff_ratio(nbs_pop_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(nbs_pop_final)
pp_check(nbs_pop_final, type = "ecdf_overlay")
pp_check(nbs_pop_final, type = "stat", stat = "mean")
pp_check(nbs_pop_final, type = "stat", stat = "min")
pp_check(nbs_pop_final, type = "stat", stat = "max")

#Marginal posterior predictive check: Pit plot
ppc_loo_pit_qq(y = nbs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(nbs_pop_final),
               lw = weights(loo1$psis_object)) 

#--------------------------------------------------------------------------------
#Extract and plot conditional effects of each predictor from best model (i.e. posterior distributions of conditional means)
#conditioning on the mean for all other predictors, yr/site effects ignored 

#Size effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
## 95% CI
ce1s_1 <- conditional_effects(nbs_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(nbs_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(nbs_pop_final , effect = "cw", re_formula = NA,
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
  geom_rug(data = nbs.dat, aes(x = cw, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data 
  labs(x = "Carapace Width (mm)", y = "") +
  theme_minimal() +
  ylim(0,215) +
  theme(axis.text=element_text(size=8)) -> sizeplot_nbs

##Julian Day
## 95% CI
ce1s_1 <- conditional_effects(nbs_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(nbs_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(nbs_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$julian
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$julian[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$julian[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$julian[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$julian[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_rug(data = nbs.dat, aes(x = julian, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Day of Year", y = "") +
  theme_minimal() +
  ylim(0,215) +
  theme(axis.text=element_text(size=8)) -> dayplot_nbs

##Snow Crab Density 
## 95% CI
ce1s_1 <- conditional_effects(nbs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(nbs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(nbs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
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
  geom_line(size = 1, color = "#084594") +
  geom_rug(data = nbs.dat, aes(x = fourth.root.cpue, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Snow Crab Density\n(Fourth root CPUE)", y = "") +
  theme_minimal() +
  ylim(0,215) -> cpueplot_nbs

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(nbs_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(nbs_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(nbs_pop_final , effect = "temperature", re_formula = NA,
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
  geom_line(size = 1, color = "#084594") +
  geom_rug(data = nbs.dat, aes(x = temperature, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = expression("Temperature " ( degree~C)), y = "Energetic Condition\n(total FA/DWT)") +
  theme_minimal() +
  ylim(0, 215) -> tempplot_nbs

#See EBS drivers script for figure compilation 
