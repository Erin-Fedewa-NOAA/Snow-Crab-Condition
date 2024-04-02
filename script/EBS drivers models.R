#Investigate drivers of condition in EBS snow crab using Bayesian multivariate models
  #a) Population-level effects of temperature and snow crab density across years
  #b) Conditional effects of temperature and density within years 

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
  #2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
  #used as response variable in all further analyses. 

# Author: Erin Fedewa
# last updated: 3/12/24

#install github version of brms with fix for truncated skew normal model runs
#if (!requireNamespace("remotes")) {
  #install.packages("remotes")
#}
#remotes::install_github("paul-buerkner/brms")

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
#look at distributions of cpue covariates 
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
         maturity != 1,
         Total_FA_Conc_WWT > 0) %>%
  mutate(year = as.factor(year),
         sex = as.factor(sex),
         region = as.factor(sample_region),
         station = as.factor(gis_station),
         temperature = as.numeric(gear_temperature),
         fourth.root.cpue = as.numeric(cpue^0.25),
         fourth.root.invert = as.numeric(total_benthic_cpue^0.25),
         julian = as.numeric(julian)) -> ebs.dat 

y_obs <- ebs.dat$Total_FA_Conc_WWT

#Assess collinearity b/w covariates 
ebs.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            fourth.root.cpue = mean(fourth.root.cpue), 
            fourth.root.invert = mean(fourth.root.invert), 
            julian = mean(julian)) -> corr.dat

cor(corr.dat[,3:6]) #All < 0.6
corrplot(cor(corr.dat[,3:6]), method = 'number') 

##########################################
#Distribution of response variable - choosing a brms family 
ebs.dat %>%
  ggplot(aes(Total_FA_Conc_WWT)) + 
  geom_density() #pretty darn left skewed 

#Test glm model to look at distribution of residuals 
test.1 <- glm(Total_FA_Conc_WWT ~ year, data=ebs.dat, family =gaussian(link=log))
plot(test.1) #overdispersion in qqplot, not gaussian! 
plot(density(resid(test.1, type='deviance'))) #very long tail, much heavier than Gaussian

# fitting a test brms model with a Gaussian likelihood - truncating response at 0, else models predict values < 0
  #ie. responses out of bounds are discarded
model_normal_trunc <- brm(Total_FA_Conc_WWT | trunc(lb = 0) ~ year, family = gaussian(link="identity"), data = ebs.dat)

#Gaussian with log response so no < 0 predictions
model_normal_log <- brm(log(Total_FA_Conc_WWT) ~ year, family = gaussian(link="identity"), data = ebs.dat)
summary(model_normal_log)

# Student t distribution (more robust to outliers)- truncating response at 0 
model_student_trunc <- brm(Total_FA_Conc_WWT | trunc(lb = 0)  ~ year, family = student(), data = ebs.dat)
#Compiling issue that I can't sort out! 

# fitting a test brms model with a skew normal likelihood (more flexible distribution)- truncating response
model_skew_trunc <- brm(Total_FA_Conc_WWT | trunc(lb = 0)  ~ year, family = skew_normal(), data = ebs.dat)

#skew normal with log response 
model_skew_log <- brm(log(Total_FA_Conc_WWT)  ~ year, family = skew_normal(), data = ebs.dat)

# fitting a test brms model with a gamma likelihood
model_gamma <- brm(Total_FA_Conc_WWT ~ year, family = "gamma", data = ebs.dat)
summary(model_gamma)

# fitting a test brms model with a lognormal likelihood
model_log <- brm(Total_FA_Conc_WWT ~ year, family = lognormal(), data = ebs.dat)
summary(model_log)

# posterior predictive checking
pp_check(model_normal_trunc, ndraws = 1e2) + pp_check(model_normal_log, ndraws = 1e2) + 
  pp_check(model_student_trunc, ndraws = 1e2) + pp_check(model_skew_trunc, ndraws = 1e2) +
  pp_check(model_skew_log, ndraws = 1e2)+ pp_check(model_gamma, ndraws = 1e2) + pp_check(model_log, ndraws = 1e2) 
#truncated skew normal captures mean and variance best, none of the models picking up apparent bimodality of data 

# posterior predictive checking - boxplots
pp_check(model_normal_trunc, type = "boxplot", ndraws = 20) + pp_check(model_normal_log, type = "boxplot", ndraws = 20) + 
  pp_check(model_student_trunc, type = "boxplot", ndraws = 20) + pp_check(model_skew_trunc, type = "boxplot", ndraws = 20) +
  pp_check(model_skew_log, type = "boxplot", ndraws = 20) + pp_check(model_gamma, type = "boxplot", ndraws = 20) + 
  pp_check(model_log, type = "boxplot", ndraws = 20) 

#let's look at the distribution of minimum values for posterior distributions vrs data
pp_check(model_normal_trunc, type = "stat", stat = "min") + pp_check(model_normal_log, type = "stat", stat = "min") + 
  pp_check(model_student_trunc, type = "stat", stat = "min") + pp_check(model_skew_trunc, type = "stat", stat = "min") +
  pp_check(model_skew_log, type = "stat", stat = "min") + pp_check(model_gamma, type = "stat", stat = "min") + 
  pp_check(model_log, type = "stat", stat = "min")

#now maximum values
pp_check(model_normal_trunc, type = "stat", stat = "max") + pp_check(model_normal_log, type = "stat", stat = "max") + 
  pp_check(model_student_trunc, type = "stat", stat = "max") + pp_check(model_skew_trunc, type = "stat", stat = "max") +
  pp_check(model_skew_log, type = "stat", stat = "max") + pp_check(model_gamma, type = "stat", stat = "max") + 
  pp_check(model_log, type = "stat", stat = "max")
#lognormal is way overshooting max 

#and means
pp_check(model_normal_trunc, type = "stat", stat = "mean") + pp_check(model_normal_log, type = "stat", stat = "mean") + 
  pp_check(model_student_trunc, type = "stat", stat = "mean") + pp_check(model_skew_trunc, type = "stat", stat = "mean") +
  pp_check(model_skew_log, type = "stat", stat = "mean") + pp_check(model_gamma, type = "stat", stat = "mean") + 
  pp_check(model_log, type = "stat", stat = "mean")
  
#Which models are providing best predictions?
model_normal_trunc <- add_criterion(model_normal_trunc, "waic")
model_student_trunc <- add_criterion(model_student_trunc, "waic")
model_skew_trunc <- add_criterion(model_skew_trunc, "waic")
model_gamma <- add_criterion(model_gamma, "waic")
model_log <- add_criterion(model_log, "waic")
  loo_compare(model_normal_trunc, model_student_trunc, model_skew_trunc, model_gamma, model_log, criterion = "waic")
#predictive accuracy highest for skew normal truncated model
  
#Let's include log(y) models, but we need Jacobian correction to make models comparable 
loo_normal_log <- loo(model_normal_log)
loo_normal_trunc <- loo(model_normal_trunc)
loo_skew_trunc <- loo(model_skew_trunc)
loo_skew_log <- loo(model_skew_log)

loo_normal_log_jacobian <- loo_normal_log
loo_normal_log_jacobian$pointwise[,1] <- loo_normal_log_jacobian$pointwise[,1] - log(ebs.dat$Total_FA_Conc_WWT)
loo_skew_log_jacobian <- loo_skew_log
loo_skew_log_jacobian$pointwise[,1] <- loo_skew_log_jacobian$pointwise[,1] - log(ebs.dat$Total_FA_Conc_WWT)

loo_compare(loo_normal_trunc, loo_skew_trunc, loo_normal_log_jacobian, loo_skew_log_jacobian)

#Overall summary: The more flexible skew normal model continues to fit the data best, and 
  #visually, a truncated skew normal model appears like the best choice
#HOWEVER, brms cannot yet compute posterior draws for truncated skew normal models yet
#Results from models with more covariates and truncated Gaussian distribution appear robust to the 
  #two distributions, and ppc plots don't look terrible, so let's go with Gaussian 
#############################################
#EBS Models: 
#Model runs not shown here, but group-level effects structure was explored. Due to the high number of stations
  #containing only 1 crab, 1|station and nested 1|region/station models had convergence issues. 
#We'll go with 1|region to at least attempt to account for the repeat sampling design, and note that
  #we're not using 1|year/region b/c this confounds covariate effects with strong annual signals like temp! 

####################################
#Goal #1: Interpret population-level effects of temperature and snow crab density across years

#MODEL 1 BASE MODEL: default, truncated Gaussian, 
  #Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)
    #Note: weakly informative priors were tested, but don't improve pp checks

mod1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      (1 | region)) 

## Show default priors
get_prior(mod1_formula, ebs.dat)

mod1 <- brm(mod1_formula,
            data = ebs.dat,
            family = gaussian,
            #prior = c(prior(student_t(3, 96, 59), class = Intercept, lb = 0),
              #prior(cauchy(0, 20),  class = sigma, lb = 0)),
                        cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod1, file = "./output/mod1.rds")
mod1 <- readRDS("./output/mod1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod1$fit)
neff_lowest(mod1$fit)
rhat_highest(mod1$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod1, ask = FALSE)
plot(conditional_effects(mod1), ask = FALSE)
mcmc_plot(mod1, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod1)) #Effective sample size: All ratios > 0.1
pp_check(mod1)

summary(mod1) #credible intervals for spline variance parameters (sds) don't include 0, let's keep smooths
bayes_R2(mod1) #R2 = 0.20
loo(mod1) -> plot(a)

###########################
#MODEL 2 BASE MODEL + INVERT: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert main effect 

mod2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k = 4) + (1 | region))  

mod2 <- brm(mod2_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))


#Save output
saveRDS(mod2, file = "./output/mod2.rds")
mod2 <- readRDS("./output/mod2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod2$fit)
neff_lowest(mod2$fit)
rhat_highest(mod2$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod2, ask = FALSE)
plot(conditional_effects(mod2), ask = FALSE)
mcmc_plot(mod2, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod2)) #Effective sample size: All ratios > 0.1
pp_check(mod2)

summary(mod2) 
bayes_R2(mod2) #R2 = 0.22
loo(mod2) -> b
plot(b)

# model comparison
loo(mod1, mod2, moment_match = TRUE) 
#Interesting...at medium to high levels of benthic prey, we see declines in 
  #energetic condition. Maybe more a proxy for competition than prey? 
#With a SE < 5, benthic inverts really doesn't add much to predictive capacity. 
  #We'll drop for parsimony until we can better resolve a benthic prey index

#######################
#MODEL 3 BASE MODEL + CRAB DENSITY: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + crab density fixed effect

mod3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                       s(fourth.root.cpue, k=4) + (1 | region))  

mod3 <- brm(mod3_formula,
            data = ebs.dat,
            family = gaussian,
            prior = c(prior(student_t(3, 96, 59), class = Intercept, lb = 0),
            prior(cauchy(0, 20),  class = sigma, lb = 0)),
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod3, file = "./output/mod3.rds")
mod3 <- readRDS("./output/mod3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod3$fit)
neff_lowest(mod3$fit)
rhat_highest(mod3$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod3, ask = FALSE)
plot(conditional_effects(mod3), ask = FALSE)
mcmc_plot(mod3, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod3)) #Effective sample size: All ratios > 0.1
pp_check(mod3)
pp_check(mod3, type = "stat_2d")

summary(mod3) 
bayes_R2(mod3) #R2 = 0.21
loo(mod3) -> c
plot(c)

# model comparison
loo(mod1, mod3, moment_match = TRUE) 
#Also pretty small difference in models (eldp_diff/SE_diff < 2, elpd_diff < 4) with an unclear
  #effect of density. We'll keep cpue in, recognizing that models are very similar 

################################
#MODEL 4 BASE MODEL + TEMP: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature fixed effect 

mod4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                       s(temperature, k=4) + (1 | region))  

mod4 <- brm(mod4_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod4, file = "./output/mod4.rds")
mod4 <- readRDS("./output/mod4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod4$fit)
neff_lowest(mod4$fit)
rhat_highest(mod4$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod4, ask = FALSE)
plot(conditional_effects(mod4), ask = FALSE)
mcmc_plot(mod4, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod4)) #Effective sample size: All ratios > 0.1
pp_check(mod4)

summary(mod4) 
bayes_R2(mod4) #R2 = 0.24
loo(mod4) -> d
plot(d)

# model comparison
loo(mod1, mod3, mod4, moment_match = TRUE)
#Temperature model has highest predictive capacity 

####################################
#MODEL 5 BASE MODEL + TEMP + DENSITY: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature and density fixed effect 

mod5_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(temperature, k=4) + s(fourth.root.cpue, k=4) +
                      (1 | region))  

mod5 <- brm(mod5_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod5, file = "./output/mod5.rds")
mod5 <- readRDS("./output/mod5.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod5$fit)
neff_lowest(mod5$fit)
rhat_highest(mod5$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod5, ask = FALSE)
plot(conditional_effects(mod5), ask = FALSE)
mcmc_plot(mod5, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod5)) #Effective sample size: All ratios > 0.1
pp_check(mod5)

summary(mod5) 
bayes_R2(mod5) #R2 = 0.25
loo(mod5) -> d
plot(d)

# model comparison
loo(mod1, mod3, mod4, mod5, moment_match = TRUE)
#So seems that full additive model has highest predictive capacity 
  #and is a substantial improvement over mod3 and mod1

###################################
#Full Model Comparison (base model + base/cpue + base/cpue/temp)

#LOO-CV
mod1 <- add_criterion(mod1, "loo")
mod3 <- add_criterion(mod3, "loo")
mod5 <- add_criterion(mod5, "loo")
loo_compare(mod1, mod3, mod5, criterion = "loo") %>% print(simplify = F)
  model.comp <- loo(mod1, mod3, mod5, moment_match = TRUE)

#and loo weights
model_weights(mod1, mod3, mod5, weights = "loo") %>% round(digits = 2)
#Again, full model is best

#Table of Rsq Values 
rbind(bayes_R2(mod1), 
      bayes_R2(mod3), 
      bayes_R2(mod5)) %>%
  as_tibble() %>%
  mutate(model = c("mod1", "mod3", "mod5"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#Model weights 
loo1 <- loo(mod1)
loo3 <- loo(mod3)
loo5 <- loo(mod5)

loo_list <- list(loo1, loo3, loo5)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
  #Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)
#Full model is consistently highest weighted model

#Save model output 
tab_model(mod1, mod3, mod5)

forms <- data.frame(formula=c(as.character(mod1_formula)[1],
                              as.character(mod3_formula)[1],
                              as.character(mod5_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/ebs_pop_model_comp.csv")

#################################
#FINAL MODEL:  Run mod5 model with 10,000 iterations and set seed for reproducibility 
ebs_pop_final <- brm(mod5_formula,
                     data = ebs.dat,
                     family = gaussian,
                     cores = 4, chains = 4, iter = 10000, warmup = 1000,
                     save_pars = save_pars(all = TRUE), seed = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(ebs_pop_final, file = "./output/ebs_pop_final.rds")
ebs_pop_final <- readRDS("./output/ebs_pop_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(ebs_pop_final$fit)
neff_lowest(ebs_pop_final$fit)
rhat_highest(ebs_pop_final$fit)
summary(ebs_pop_final)
bayes_R2(ebs_pop_final) #r2 = .24
loo(ebs_pop_final)

#Diagnostic Plots
plot(ebs_pop_final, ask = FALSE)
plot(conditional_effects(ebs_pop_final), ask = FALSE)
mcmc_plot(ebs_pop_final, prob = 0.95)
mcmc_neff(neff_ratio(ebs_pop_final)) #Effective sample size: All ratios > 0.1
hypothesis(ebs_pop_final, "stemperature_1" < 0)

#Posterior Predictive Check Plots:
pp_check(ebs_pop_final)
pp_check(ebs_pop_final, type = "ecdf_overlay")
pp_check(ebs_pop_final, type = "stat", stat = "mean")
pp_check(ebs_pop_final, type = "stat", stat = "min")
pp_check(ebs_pop_final, type = "stat", stat = "max")
  
#Pit plots
pit <- function(y, yrep) {
  n_draws <- nrow(yrep)
  pit <- sapply(1:length(y),
                \(n) {
                  mean(y[n] > yrep[, n]) +
                    # randomized PIT for discrete y (Czado, C., Gneiting, T.,
                    # Held, L.: Predictive model assessment for count
                    # data. Biometrics 65(4), 1254–1261 (2009).)
                    sample(sum(y[n] == yrep[, n]), 1) / n_draws
                })
  pmax(pmin(pit, 1), 0)
}

ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_pop_final))) #no overdisersion, looks good
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_pop_final))

################################
#Extract and plot conditional effects of each predictor from best model (i.e. posterior distributions of conditional means)
  #conditioning on the mean for all other predictors, yr/site effects ignored 

#tidybayes method: massive dataset being passed to functions crashing R....skip to line 522

#Plot posterior distributions of conditional means 
ebs.dat %>%
  #generate grid with temperature predictions
  data_grid(temperature = seq_range(temperature, n=100)) %>%
  #add draws from posterior distributions of conditional means
  add_epred_draws(ebs_pop_final, re_formula = NA) -> dat.epred #no group level effects

#temperature
dat.epred %>%
  ggplot(aes(x = temperature, y = Total_FA_Conc_WWT)) +
   stat_lineribbon(aes(y = .epred)) +
    geom_point(data = ebs.dat) 
    
#Plot posterior predictions
ebs.dat %>%
  data_grid(temperature, cw, julian, fourth.root.cpue) %>%
  add_predicted_draws(ebs_pop_final, re_formula = NA) -> dat.pospred

#temperature
dat.pospred %>%
  ggplot(aes(x = temperature, y = Total_FA_Conc_WWT)) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80), alpha = 1/4) +
  geom_point(data = ebs.dat) 

 #Size effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
## 95% CI
ce1s_1 <- conditional_effects(ebs_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(ebs_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_pop_final , effect = "cw", re_formula = NA,
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
  geom_rug(data = ebs.dat, aes(x = cw, y = Total_FA_Conc_WWT), 
           colour = "grey80", size = .75, sides="b") + #raw data) 
  labs(x = "Carapace width (mm)", y = "Energetic Condition (Total FA/WWT)") +
  theme_bw() +
  ylim(0,225) -> sizeplot

##Julian Day
## 95% CI
ce1s_1 <- conditional_effects(ebs_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(ebs_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_pop_final , effect = "julian", re_formula = NA,
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
  geom_rug(data = ebs.dat, aes(x = julian, y = Total_FA_Conc_WWT), 
           colour = "grey80", size = .75, sides="b") + #raw data) 
  labs(x = "Day of Year", y = "") +
  theme_bw() +
  ylim(0,225) -> dayplot

##Snow Crab Density 
## 95% CI
ce1s_1 <- conditional_effects(ebs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(ebs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_pop_final , effect = "fourth.root.cpue", re_formula = NA,
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
  geom_rug(data = ebs.dat, aes(x = fourth.root.cpue, y = Total_FA_Conc_WWT), 
           colour = "grey80", size = .75, sides="b") + #raw data
  labs(x = "Snow Crab Density (4th root CPUE)", y = "Energetic Condition (Total FA/WWT)") +
  theme_bw() +
  ylim(0, 225) -> cpueplot

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(ebs_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(ebs_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_pop_final , effect = "temperature", re_formula = NA,
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
  geom_rug(data = ebs.dat, aes(x = temperature, y = Total_FA_Conc_WWT), colour = "grey80", 
           size = .75, sides="b") + #raw data
  labs(x = "Temperature (C)", y = "") +
  theme_bw() +
  ylim(0, 225) -> tempplot

#Combine plots for Fig 4 of MS
(sizeplot + dayplot) / (cpueplot + tempplot) + 
  plot_annotation(tag_levels = 'a', title = "Eastern Bering Sea Snow Crab",
        theme = theme(plot.title = element_text(hjust = 0.5)))
ggsave("./figs/ebs_pop_Fig4.png")

#####################################################
#Marginal Effects: instantaneous slope of one explanatory value with all 
#other values held constant

#Marginal effect at the mean: julian day slope
ebs_pop_final %>%
  emtrends(~ julian, 
           var = "julian", 
           regrid = "response")
#on average, a one-day increase in Julian day is associated with a 1.1% increase in 
#the probability of infection

#Marginal effect at various levels of julian day  
ebs_pop_final %>% 
  emtrends(~ julian, var = "julian",
           at = list(julian = 
                       seq(min(ebs.dat$julian), 
                           max(ebs.dat$julian), 1)),
           re_formula = NA) %>%
  as_tibble() %>%
  #and plot 
  ggplot(aes(x = julian, y = julian.trend)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.1) +
  geom_line(size = 1) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Julian Day", y = "Marginal effect of julian day on probability of infection") +
  theme_bw() 

#Marginal effect at the mean: cw 
ebs_pop_final %>%
  emtrends(~ cw, 
           var = "cw", 
           regrid = "response", re_formula = NA)
#a 1mm increase in Julian day is associated with a 1.1% increase in 
#the probability of infection

#Marginal effect at various size crab 
ebs_pop_final %>% 
  emtrends(~ cw, var = "cw",
           at = list(cw = 
                       seq(min(ebs.dat$cw), 
                           max(ebs.dat$cw), 1)),
           re_formula = NA) %>%
  as_tibble() %>%
  #and plot 
  ggplot(aes(x = cw, y = cw.trend)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.1) +
  geom_line(size = 1) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Carapace width", y = "Marginal effect of size on probability of infection") +
  theme_bw()  

######################################################
#Generating posterior predictions for final model 

#global size mean-ignoring year/site specific deviations 
grand_mean <- ebs_pop_final %>% 
  #create dataset across a range of observed sizes sampled
  epred_draws(newdata = expand_grid(size = range(ebs.dat$size),
                                    temperature = mean(ebs.dat$temperature), 
                                    julian = mean(ebs.dat$julian)), 
              re_formula = NA) #ignoring random effects 
#plot
ggplot(grand_mean, aes(x = size, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Carapace width", y = "Probability of infection",
       fill = "Credible interval") +
  theme_bw() +
  theme(legend.position = "bottom")

#average marginal effect of size: i.e. finding the slope at different sizes 
grand_mean_ame <- ebs_pop_final %>% 
  emtrends(~ size,
           var = "size",
           at = list(julian = mean(ebs.dat$julian),
                     temperature=mean(ebs.dat$temperature),
                     size = c(30, 60, 90)),
           epred = TRUE, re_formula = NA) %>% 
  #get predicted values from posterior draws 
  gather_emmeans_draws()

ggplot(grand_mean_ame, aes(x = .value, fill = factor(size))) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of an increase in crab size",
       y = "Density", fill = "Size") +
  theme_bw() 
#Sampling a 30mm crab is associated with a ~1% increase in prob of infection- 
#smaller the size, larger the marginal effect 

#Average overall slope at mean size 
ebs_pop_final %>% 
  emtrends(~ 1,
           var = "size",
           epred = TRUE, re_formula = NA) 

#####

#Year-specific posterior predictions across size 
all_years <- ebs_pop_final %>% 
  epred_draws(newdata = expand_grid(size = range(ebs.dat$size),
                                    temperature = mean(ebs.dat$temperature), 
                                    julian = mean(ebs.dat$julian), 
                                    year = levels(ebs.dat$year)), 
              re_formula = ~ (1 | year)) #only predict using yr effects, not site too 

ggplot(all_years, aes(x = size, y = .epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Carapace width", y = "Probability of Infection",
       fill = "Credible interval") +
  facet_wrap(vars(year)) +
  theme_bw() +
  theme(legend.position = "bottom")

#average marginal effect by year
all_years_ame <- ebs_pop_final %>% 
  emtrends(~ size + year,
           var = "size",
           at = list(year = levels(ebs.dat$year)),
           epred = TRUE, re_formula = ~ (1 | year)) %>% 
  gather_emmeans_draws()

ggplot(all_years_ame,aes(x = .value)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect of a\1-point increase in crab size",
       y = "Density") +
  facet_wrap(~year) +
  theme_bw()

#post and interval summaries of draws from size effect 
all_years_ame %>% median_hdi()
#Very little variation in size effect across years 

####################################
#Goal #2: Interpret conditional effects of temperature and snow crab density within each year
  #Testing best model from Goal #1 with year interactions

#MODEL 6 BASE MODEL + TEMP*YR + DENSITY*YR: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature*year and 
  #density*year fixed effect 

mod6_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(temperature, k=4, by=year) + s(fourth.root.cpue, k=4, by=year) +
                      (1 | region))  

mod6 <- brm(mod6_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(mod6, file = "./output/mod6.rds")
mod6 <- readRDS("./output/mod6.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod6$fit)
neff_lowest(mod6$fit)
rhat_highest(mod6$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(mod6, ask = FALSE)
plot(conditional_effects(mod6), ask = FALSE)
mcmc_plot(mod6, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod6)) #Effective sample size: All ratios > 0.1
pp_check(mod6)

summary(mod6) 
bayes_R2(mod6) #R2 = 0.45
loo(mod6) -> d
plot(d)

# model comparison
loo(mod5, mod6, moment_match = TRUE)

# model comparison
mod1 <- add_criterion(mod1, "loo")
mod3 <- add_criterion(mod3, "loo")
mod4 <- add_criterion(mod4, "loo")
mod5 <- add_criterion(mod5, "loo")
mod6 <- add_criterion(mod6, "loo")
loo_compare(mod1, mod3, mod4, mod5, mod6, criterion = "loo") %>% print(simplify = F)

#and loo weights
model_weights(mod1, mod3, mod4, mod5, mod6, weights = "loo") %>% round(digits = 2)
#So seems that full model with interactions has highest predictive capacity

#################################
#FINAL MODEL:  Run mod6 model with 10,000 iterations and set seed for reproducibility 
ebs_yrixn_final <- brm(mod6_formula,
                     data = ebs.dat,
                     family = gaussian,
                     cores = 4, chains = 4, iter = 10000, warmup = 1000,
                     save_pars = save_pars(all = TRUE), seed = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(ebs_yrixn_final, file = "./output/ebs_yrixn_final.rds")
ebs_yrixn_final <- readRDS("./output/ebs_yrixn_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(ebs_yrixn_final$fit)
neff_lowest(ebs_yrixn_final$fit)
rhat_highest(ebs_yrixn_final$fit)
summary(ebs_yrixn_final)
bayes_R2(ebs_yrixn_final) 
loo(ebs_yrixn_final)

#Diagnostic Plots
plot(ebs_yrixn_final, ask = FALSE)
plot(conditional_effects(ebs_yrixn_final), ask = FALSE)
mcmc_plot(ebs_yrixn_final, prob = 0.95)
mcmc_neff(neff_ratio(ebs_yrixn_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(ebs_yrixn_final) #this doesn't look great...
pp_check(ebs_yrixn_final, type = "ecdf_overlay")
pp_check(ebs_yrixn_final, type = "stat", stat = "mean")
pp_check(ebs_yrixn_final, type = "stat", stat = "min")
pp_check(ebs_yrixn_final, type = "stat", stat = "max")

#Pit plots
pit <- function(y, yrep) {
  n_draws <- nrow(yrep)
  pit <- sapply(1:length(y),
                \(n) {
                  mean(y[n] > yrep[, n]) +
                    # randomized PIT for discrete y (Czado, C., Gneiting, T.,
                    # Held, L.: Predictive model assessment for count
                    # data. Biometrics 65(4), 1254–1261 (2009).)
                    sample(sum(y[n] == yrep[, n]), 1) / n_draws
                })
  pmax(pmin(pit, 1), 0)
}

ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_yrixn_final))) #slight overdispersion
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_yrixn_final))

################################
#Extract and plot conditional effects of yr*cpue and yr*temperature interaction





#Priors? interaction model doesn't look great...
#ditch table/model selection altogether? so phrase not as trying to determine drivers,
  #but evaluate the effects of temp and density, phrasing intro as collapse/heat wave 
#Sort out conditional vrs marginal - using the correct for Figs?
#run final interaction model
#figure out how to plot interactions for cpue/temp
#seperate figure or combine?
#run same models for NBS


#overall temp and density effect on conditon (avg effect across years)
#And then, look at year interaction with best fit model to see if these 
  #variables differ between pre and post collapse yrs (within yrs)
#The strength of association and direction b/w temperature differs in heat wave yr
#so strong temp effect, but yr still explains much of variation..likely b/c temp
  #is a proxy for poor ecosystem conditions during collapse/heat wave

#Read chp 8
#read visualization paper Fig 6/9/10
#use priorsense package to test prior and likelihood sensitivity
#
#Wed: 
#Need to set seed and up iteration for final run
#any other model diagnostics/comparisons
#Extract draws and make plots - how to interpret interactions vrs overall effect?

pit <- function(y, yrep) {
  n_draws <- nrow(yrep)
  pit <- sapply(1:length(y),
                \(n) {
                  mean(y[n] > yrep[, n]) +
                    # randomized PIT for discrete y (Czado, C., Gneiting, T.,
                    # Held, L.: Predictive model assessment for count
                    # data. Biometrics 65(4), 1254–1261 (2009).)
                    sample(sum(y[n] == yrep[, n]), 1) / n_draws
                })
  pmax(pmin(pit, 1), 0)
}

ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(mod2))) #overdispersion
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(mod2))


#try this - posterior predictions via tidybayes
ebs.dat %>%
  group_by(year) %>%
  data_grid(temperature = seq_range(temperature, n = 101)) %>%
  add_predicted_draws(mod2, re_formula = NA, category="temperature") %>% #can use add_epred_draws() instead for posterior means
  ggplot(aes(x = temperature, y = Total_FA_Conc_WWT, color = ordered(year), fill = ordered(year))) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50), alpha = 1/4) +
  geom_point(data = ebs.dat, colour = "darkseagreen4", size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Dark2")

#or faceted
#facet_grid(.~year, space="free_x", scales="free_x")
ame_fancy_zi_quota <- mod2 %>%
  avg_comparisons(variables = "temperature") %>% 
  posterior_draws()

ggplot(ame_fancy_zi_quota, aes(x = draw)) +
  stat_halfeye(.width = c(0.8, 0.95), point_interval = "median_hdi",
               fill = "#bc3032") +
  labs(x = "Average marginal effect of having a gender-based\nquota on the proportion of women in parliament", y = NULL,
       caption = "80% and 95% credible intervals shown in black") +
  theme_clean()

r_fancy <- ame_fancy_zi_quota %>% median_hdi(draw) #after accounting for other covariates, FA conc
#of snow crab 

mcmc_areas(as.matrix(mod2), regex_pars = "temperature")



#Add a, b and c to figs
patchwork + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))



