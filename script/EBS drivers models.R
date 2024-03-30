#Investigate drivers of condition in EBS snow crab using Bayesian multivariate models

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
plot(conditional_smooths(mod1), ask = FALSE)
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
plot(conditional_smooths(mod2), ask = FALSE)
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
plot(conditional_smooths(mod3), ask = FALSE)
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
#Also pretty small difference in models (SE < 5) with an unclear
  #effect of density 

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
plot(conditional_smooths(mod4), ask = FALSE)
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
plot(conditional_smooths(mod5), ask = FALSE)
mcmc_plot(mod5, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod5)) #Effective sample size: All ratios > 0.1
pp_check(mod5)

summary(mod5) 
bayes_R2(mod5) #R2 = 0.24
loo(mod5) -> d
plot(d)

# model comparison
loo(mod1, mod3, mod4, mod5, moment_match = TRUE)
#So seems that full model has highest predictive capacity 

####################################
#MODEL 6 BASE MODEL + TEMP * DENSITY: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature*density interaction 

mod6_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                    s(fourth.root.cpue, temperature, k=4) + (1 | region))  

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
plot(conditional_smooths(mod6), ask = FALSE)
mcmc_plot(mod6, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod6)) #Effective sample size: All ratios > 0.1
pp_check(mod6)

summary(mod6) 
bayes_R2(mod6) #R2 = 0.24
loo(mod6) -> d
plot(d)

# model comparison
loo(mod1, mod3, mod4, mod5, mod6, moment_match = TRUE)
#Hmmm no real advantage of an interaction?

#Follow ups: How do these models compare to ones with year interaction? 
  #and is this actually a different question? I think so
#Look into model selection for interactions, brms - such low Rsq with no year effect 
  #so strong temp effect, but yr still explains much of variation 







###################################
#Full Model Comparison

#LOO-CV
model.comp <- loo(mod1, mod2, mod3, mod4, moment_match = TRUE)
model.comp
#Again, full model is best

#Table of Rsq Values 
rbind(bayes_R2(mod1), 
      bayes_R2(mod2), 
      bayes_R2(mod3), 
      bayes_R2(mod4)) %>%
  as_tibble() %>%
  mutate(model = c("mod1", "mod2", "mod3", "mod4"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#Model weights 
#PSIS-LOO
loo1 <- loo(mod1)
loo2 <- loo(mod2)
loo3 <- loo(mod3)
loo4 <- loo(mod4)

loo_list <- list(loo1, loo2, loo3, loo4)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
#Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)
#Tanner3 consistently highest weighted model

#Save model output 
tab_model(mod1, mod2, mod3, mod4)

forms <- data.frame(formula=c(as.character(mod1_formula)[1],
                              as.character(mod2_formula)[1],
                              as.character(mod3_formula)[1],
                              as.character(mod4_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/ebs_model_comp.csv")

#################################
#Final model using tanner3, though all models tested are very similar

#Final Model:  Run tanner3 model with 10,000 iterations and set seed for reproducibility 
tannerfinal <- brm(tanner3_formula,
                   data = tanner.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, iter = 10000,
                   save_pars = save_pars(all = TRUE), seed = 3, 
                   control = list(adapt_delta = 0.9999, max_treedepth = 14))

#Save model output 
saveRDS(tannerfinal, file = "./output/tannerfinal.rds")
tannerfinal <- readRDS("./output/tannerfinal.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tannerfinal$fit)
neff_lowest(tannerfinal$fit)
rhat_highest(tannerfinal$fit)
summary(tannerfinal)
bayes_R2(tannerfinal)

#Diagnostic Plots
plot(tannerfinal, ask = FALSE)
plot(conditional_smooths(tannerfinal), ask = FALSE)
mcmc_plot(tannerfinal, prob = 0.95)
mcmc_plot(tannerfinal, transformations = "inv_logit_scaled")
mcmc_rhat(rhat(tannerfinal)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tannerfinal, pars = c("b_Intercept", "bs_ssize_1", "bs_stemperature_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tannerfinal)) #Effective sample size: All ratios > 0.1
marginal_effects(tannerfinal, surface = TRUE) #visualize effects of predictors on the expected response
marginal_smooths(tannerfinal) #
hypothesis(tannerfinal, "ssize_1 < 0")

#Posterior Predictive Check: Mean and skewness summary statistics 
color_scheme_set("red")
pmean1 <- posterior_stat_plot(y_obs, tannerfinal) + 
  theme(legend.text = element_text(size=8), 
        legend.title = element_text(size=8)) +
  labs(x="Mean", title="Mean")

color_scheme_set("gray")
pskew1 <- posterior_stat_plot(y_obs,tannerfinal, statistic = "skew") +
  theme(legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  labs(x = "Fisher-Pearson Skewness Coeff", title="Skew")

pmean1 + pskew1

#PPC: Classify posterior probabilities and compare to observed 
preds <- posterior_epred(tannerfinal)
pred <- colMeans(preds) #averaging across draws 
pr <- as.integer(pred >= 0.5) #Classify probabilities >0.5 as presence of disease 
mean(xor(pr, as.integer(y_obs == 0))) # posterior classification accuracy looks good

# Compute AUC for predicting prevalence with the model
y_obs <- tanner.dat$pcr
preds <- posterior_epred(tannerfinal)
auc <- apply(preds, 1, function(x) {
  roc <- roc(y_obs, x, quiet = TRUE)
  auc(roc)
})
hist(auc) #Looks like our model discriminates fairly well 


################################
#Extract and plot conditional effects of each predictor from best model
#conditioning on the mean for all other predictors, yr/site effects ignored 

#Size effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "size", re_formula = NA,
                              probs = c(0.1, 0.9))

dat_ce <- ce1s_1$size
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$size[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$size[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$size[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$size[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), fill = "#F7FBFF") +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), fill = "#DEEBF7") +
  geom_ribbon(aes(ymin = lower_80, ymax = upper_80), fill = "#C6DBEF") + 
  geom_line(size = 1, color = "black") +
  geom_point(data = tanner.dat, aes(x = size, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Carapace width (mm)", y = "Probability of infection") +
  theme_bw() -> sizeplot

##Julian Day
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "julian", re_formula = NA,
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
  geom_point(data = tanner.dat, aes(x = julian, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Julian Day", y = "") +
  theme_bw() -> dayplot

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(tannerfinal , effect = "temperature", re_formula = NA,
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
  geom_point(data = tanner.dat, aes(x = temperature, y = pcr), colour = "grey80", shape= 73, size = 2) + #raw data
  labs(x = "Temperature (C)", y = "") +
  theme_bw() -> tempplot

#Combine plots for Fig 6 of MS
sizeplot + dayplot + tempplot + plot_annotation(tag_levels = 'a')
ggsave("./figs/tannerFig6.png")

#####################################################
#Marginal Effects: instantaneous slope of one explanatory value with all 
#other values held constant

#Marginal effect at the mean: julian day slope
tannerfinal %>%
  emtrends(~ julian, 
           var = "julian", 
           regrid = "response")
#on average, a one-day increase in Julian day is associated with a 1.1% increase in 
#the probability of infection

#Marginal effect at various levels of julian day  
tannerfinal %>% 
  emtrends(~ julian, var = "julian",
           at = list(julian = 
                       seq(min(tanner.dat$julian), 
                           max(tanner.dat$julian), 1)),
           re_formula = NA) %>%
  as_tibble() %>%
  #and plot 
  ggplot(aes(x = julian, y = julian.trend)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.1) +
  geom_line(size = 1) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Julian Day", y = "Marginal effect of julian day on probability of infection") +
  theme_bw() 

#Marginal effect at the mean: size 
tannerfinal %>%
  emtrends(~ size, 
           var = "size", 
           regrid = "response", re_formula = NA)
#a 1mm increase in Julian day is associated with a 1.1% increase in 
#the probability of infection

#Marginal effect at various size crab 
tannerfinal %>% 
  emtrends(~ size, var = "size",
           at = list(size = 
                       seq(min(tanner.dat$size), 
                           max(tanner.dat$size), 1)),
           re_formula = NA) %>%
  as_tibble() %>%
  #and plot 
  ggplot(aes(x = size, y = size.trend)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), alpha = 0.1) +
  geom_line(size = 1) +
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Carapace width", y = "Marginal effect of size on probability of infection") +
  theme_bw()  

######################################################
#Generating posterior predictions for final model 

#global size mean-ignoring year/site specific deviations 
grand_mean <- tannerfinal %>% 
  #create dataset across a range of observed sizes sampled
  epred_draws(newdata = expand_grid(size = range(tanner.dat$size),
                                    temperature = mean(tanner.dat$temperature), 
                                    julian = mean(tanner.dat$julian)), 
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
grand_mean_ame <- tannerfinal %>% 
  emtrends(~ size,
           var = "size",
           at = list(julian = mean(tanner.dat$julian),
                     temperature=mean(tanner.dat$temperature),
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
tannerfinal %>% 
  emtrends(~ 1,
           var = "size",
           epred = TRUE, re_formula = NA) 

#####

#Year-specific posterior predictions across size 
all_years <- tannerfinal %>% 
  epred_draws(newdata = expand_grid(size = range(tanner.dat$size),
                                    temperature = mean(tanner.dat$temperature), 
                                    julian = mean(tanner.dat$julian), 
                                    year = levels(tanner.dat$year)), 
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
all_years_ame <- tannerfinal %>% 
  emtrends(~ size + year,
           var = "size",
           at = list(year = levels(tanner.dat$year)),
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

############################################################
#Lastly, to compare prob of infection amoung years, lets use best model
#with year as a fixed effect 

## fit Tanner model
tanner_year_formula <-  bf(pcr ~ s(size, k = 4) + s(julian, k = 4) + s(temperature, k = 4) + year + (1 | index))

tanner_year <- brm(tanner_year_formula,
                   data = tanner.dat,
                   family = bernoulli(link = "logit"),
                   cores = 4, chains = 4, iter = 2500,
                   save_pars = save_pars(all = TRUE),
                   control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(tanner_year, file = "./output/tanner_year.rds")
tanner_year <- readRDS("./output/tanner_year.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(tanner_year$fit)
neff_lowest(tanner_year$fit)
rhat_highest(tanner_year$fit)
summary(tanner_year)
bayes_R2(tanner_year)

#Diagnostic Plots
plot(tanner_year, ask = FALSE)
plot(conditional_smooths(tanner_year), ask = FALSE)
mcmc_plot(tanner_year, prob = 0.95)
mcmc_plot(tanner_year, transformations = "inv_logit_scaled")
mcmc_rhat(rhat(tanner_year)) #Potential scale reduction: All rhats < 1.1
mcmc_acf(tanner_year, pars = c("b_Intercept", "bs_ssize_1", "bs_stemperature_1"), lags = 10) #Autocorrelation of selected parameters
mcmc_neff(neff_ratio(tanner_year)) #Effective sample size: All ratios > 0.1
marginal_effects(tanner_year, surface = TRUE) #visualize effects of predictors on the expected response
marginal_smooths(tanner_year) #
hypothesis(tanner_year, "year2016 < 0")

#Conditional Effect 
conditional_effects(tanner_year, effect = "year")

ce1s_1 <- conditional_effects(tanner_year, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Tanner crab") -> year_tanner

#Average marginal effect of year 
years_ame <- tanner_year %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NA) %>% 
  gather_emmeans_draws()

ggplot(years_ame,aes(x = .value, fill=year)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect",
       y = "Density") +
  theme_bw()

#Combine tanner/snow effects (run lines 862-902 in analyze_opilio.R first)
dodge <- position_dodge(width=0.5) #to offset datapoints on plot 

new_colors <- c("#238b45","#2171b5")

year_tanner %>%
  full_join(year_snow) %>%
  #Combined conditional effect plot 
  ggplot() +
  geom_point(aes(year, estimate__, color=factor(species, 
                                                levels = c("Tanner crab", "Snow crab"))), size=3,
             position=dodge) +
  geom_errorbar(aes(year, ymin=lower__, ymax=upper__, color=factor(species, 
                                                                   levels = c("Tanner crab", "Snow crab"))), width=0.3, 
                size=0.5, position=dodge) +
  ylab("Probability of infection") + xlab("") +
  scale_colour_manual(values = new_colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.title= element_blank())
ggsave("./figs/annual_brm.png", dpi=300)













###Mess below this###############

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

#Add a, b and c to figs
patchwork + 
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 8))



