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
library(performance)
library(marginaleffects)
library(emmeans)
library(corrplot)
library(factoextra)
library(patchwork)
library(modelr)
library(broom.mixed)
library(pROC)
library(ggpubr)
library(interactions)
library(priorsense)
library(ggthemes)
library(tidybayes)
library(RColorBrewer)
library(bayestestR)
library(knitr)
library(loo)
library(sjPlot)
source("./script/stan_utils.R")

#load data
condition_master <- read.csv("./data/total_FA_master.csv")

#colors for plotting
my_colors <- c("#D55E00","#9ECAE1", "#4292C6", "#084594")
my_colors2 <- c("#2A788EFF", "#E7B800", "#440154FF", "#22A884FF")

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
#We'll go with 1|region to account for the repeat sampling design, and note that
  #we're not using 1|year/region b/c this confounds covariate effects with strong annual signals like temp! 

####################################
#Goal #1: Interpret population-level effects of temperature and snow crab density across years

#MODEL 1 BASE MODEL: default priors, truncated Gaussian, 
  #Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)
   
mod1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      (1 | region)) 

## Show default priors
get_prior(mod5_formula, ebs.dat)

#informing priors- scaled based on the mean and standard deviation of x and y in the population
sd_x <- sd(ebs.dat$temperature)
sd_x2 <- sd(ebs.dat$fourth.root.cpue)
sd_x3 <- sd(ebs.dat$julian)
sd_x4 <- sd(ebs.dat$cw)
sd_y <- sd(ebs.dat$Total_FA_Conc_WWT) #59 
mean_y <- mean(ebs.dat$Total_FA_Conc_WWT) #96 
sd_temp <- 2.5*sd_y/sd_x #96
sd_cpue <- 2.5*sd_y/sd_x2 #29
sd_julian <- 2.5*sd_y/sd_x3 #15
sd_cw <- 2.5*sd_y/sd_x4 #8

#Setting weak priors with heavy tails - this is for full model 5, with all covariates  
priors_set <- c(set_prior("student_t(3, 0, 96)", class = "b", coef = "stemperature_1"),
                set_prior("student_t(3, 0, 29)", class = "b", coef = "sfourth.root.cpue_1"),
                set_prior("student_t(3, 0, 15)", class = "b", coef = "sjulian_1"),
                set_prior("student_t(3, 0, 8)", class = "b", coef = "scw_1"),
                set_prior("student_t(3, 96, 59)", class = "Intercept", lb =0),
                set_prior("student_t(3, 0, 49.1)", class = "sds"), #set from get_priors() on model run with default priors 
                set_prior("student_t(3, 0, 49.1)", class = "sigma"))

#Weakly informative priors were tested, but don't improve pp checks or sensitivities,
  #we'll stick with default priors

mod1 <- brm(mod1_formula,
            data = ebs.dat,
            family = gaussian,
            #prior = priors_set,
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
bayes_R2(mod1) #R2 = 0.19
loo(mod1) -> plot(a)

###########################
#MODEL 2 BASE MODEL + INVERT: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert main effect 

mod2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(fourth.root.invert, k = 3) + (1 | region))  

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
bayes_R2(mod2) #R2 = 0.21
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

mod3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                       s(fourth.root.cpue, k=3) + (1 | region))  

mod3 <- brm(mod3_formula,
            data = ebs.dat,
            family = gaussian,
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
r2_bayes(ebs_pop_final) #conditional r2 takes both fixed and random effects into account 
loo1 <- loo(ebs_pop_final, save_psis = TRUE)
plot(loo1)

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
  
#Marginal posterior predictive checks: Pit plots
ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_pop_final))) #no overdisersion, looks good
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_pop_final))

ppc_loo_pit_qq(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_pop_final),
  lw = weights(loo1$psis_object)) #looks fairly uniform, no clear model mispecifications

#Determining the sensitivity of the posterior to perturbations of the prior and likelihood (priorsense package)
  #prior sensitivity > .05 = prior-data conflict, and likelihood sensitivity < .05 = noninformative likelihood
powerscale_sensitivity(ebs_pop_final) %>% print(n = 50)
  #prior-data conflict = mismatch between prior and observed data, priors not on appropriate scale for predictors 

#now visually inspect
pss <- powerscale_sequence(ebs_pop_final) #estimate posterior draws based on power-scaling
#estimate posterior draws based on power-scaling
powerscale_plot_ecdf(pss, variables = c("sigma", "bs_sjulian_1", "bs_sfourth.root.cpue_1", "sds_scw_1"))

#Detection of prior-data conflicts suggests a model modification using heavy tailed priors -
  #this was tried (see model 1 script) and pp_checks did not improve 

################################
#Extract and plot conditional effects of each predictor from best model (i.e. posterior distributions of conditional means)
  #conditioning on the mean for all other predictors, yr/site effects ignored 

#tidybayes method: massive dataset being passed to add_epred_draws() crashing R....skip to line 522
  #Matthew Kay advice: "if you are just creating the huge long format data frame as an intermediate 
  #step (e.g. you are summarizing it down later), one way to solve this is to split up the input 
  #prediction grid into chunks, and pass each chunk to add_epred_draws and do the summarization, 
  #then combine the summaries i.e. split the output of data_grid(), summarize output 1 with epred, 
  #then do the same with output 2 - ... Do this via a loop, or map()"

#Plot posterior distributions of conditional means 
#ebs.dat %>%
  #generate grid with temperature predictions
 # data_grid(temperature = seq_range(temperature, n=100)) %>%
  #add draws from posterior distributions of conditional means
 # add_epred_draws(ebs_pop_final, re_formula = NA) -> dat.epred #no group level effects

#temperature
#dat.epred %>%
 # ggplot(aes(x = temperature, y = Total_FA_Conc_WWT)) +
  # stat_lineribbon(aes(y = .epred)) +
   # geom_point(data = ebs.dat) 
    
#Plot posterior predictions
#ebs.dat %>%
  #data_grid(temperature, cw, julian, fourth.root.cpue) %>%
  #add_predicted_draws(ebs_pop_final, re_formula = NA) -> dat.pospred

#temperature
#dat.pospred %>%
 # ggplot(aes(x = temperature, y = Total_FA_Conc_WWT)) +
  #stat_lineribbon(aes(y = .prediction), .width = c(.95, .80), alpha = 1/4) +
  g#eom_point(data = ebs.dat) 

####
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
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data 
  labs(x = "Carapace Width", y = "") +
  theme_minimal() +
  ylim(0,215) -> sizeplot

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
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Day of Year", y = "") +
  theme_minimal() +
  ylim(0,215) -> dayplot

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
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Snow Crab Density", y = "") +
  theme_minimal() +
  ylim(0,215) -> cpueplot

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
  geom_rug(data = ebs.dat, aes(x = temperature, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Bottom Temperature", y = "") +
  theme_minimal() +
  ylim(0, 215) -> tempplot

#Combine plots for Fig 4 of MS
(sizeplot + dayplot) / (cpueplot + tempplot) + 
  plot_annotation(tag_levels = 'a', title = "Eastern Bering Sea Snow Crab",
        theme = theme(plot.title = element_text(hjust = 0.5))) -> plot

#workaround to get a single shared y axis label with patchwork 
wrap_elements(plot) +
  labs(tag = "Energetic Condition (Total FA/WWT)") +
  theme(plot.tag = element_text(size = rel(1), angle = 90),
    plot.tag.position = "left") -> final_plot
ggsave("./figures/ebs_pop_Fig4.png", plot=final_plot)

#####################################################
#Other miscellaneous snippets of code 

#probability of direction of effect 
temp_draws <- ebs_pop_final %>%
  spread_draws(bs_stemperature_1)

ggplot(temp_draws, aes(x = bs_stemperature_1)) +
  stat_halfeye() 

  #Find the proportion of posterior draws that are bigger than 0
temp_draws %>% 
  summarize(prop_greater_0 = sum(bs_stemperature_1 < 0) / n())
#99% chance a temperature effect is negative  

#Conditional predictions/effect = average region - re_formula = NA: effect of x in an average region
#random offsets of region set to 0, so ignoring them 
#Marginal predictions/effect = regions on average - re_formula = NULL: average effect of x across all regions

#average marginal effect of temperature while holding other pop-level effects at mean
pred_temp <- ebs_pop_final %>% 
  epred_draws(newdata = expand_grid(temperature = seq_range(ebs.dat$temperature, n=100), 
                                    cw = mean(ebs.dat$cw), 
                                    julian = mean(ebs.dat$julian), 
                                    fourth.root.cpue = mean(ebs.dat$fourth.root.cpue), 
                                    region = levels(ebs.dat$region)), 
              re_formula = NULL) #random effects included 

ggplot(pred_temp, aes(x = temperature, y = .epred)) +
  stat_lineribbon() + 
  scale_fill_brewer(palette = "Reds") +
  labs(x = "Temperature", y = "Predicted energetic condition",
       fill = "Credible interval") +
  theme_clean() +
  theme(legend.position = "bottom")

#Average marginal effect of temperature at different values 
ame_temp <- ebs_pop_final %>% 
  emtrends(~ temperature,
           var = "temperature",
           at = list(cw = mean(ebs.dat$cw), 
                     julian = mean(ebs.dat$julian), 
                     fourth.root.cpue = mean(ebs.dat$fourth.root.cpue),
                     temperature = c(-1, 1, 3, 5)),
           epred = TRUE) %>% 
  gather_emmeans_draws() 

ggplot(ame_temp, aes(x = .value, fill = factor(temperature))) +
  stat_halfeye(slab_alpha = 0.75) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Average marginal effect of a 1C increase in bottom temperature",
       y = "Density", fill = "Temperature (C)",
       caption = "80% and 95% credible intervals shown in black") +
  theme_clean() + 
  theme(legend.position = "bottom")
#A 2C increase in temp results in a 15% decline in condition
  #I don't think this is correct?

#And now let's play around with ROPEs: 
  #https://easystats.github.io/bayestestR/articles/region_of_practical_equivalence.html#how-to-define-the-rope-range-

perc_in_rope <- rope(ebs_pop_final, ci = 1) 
#output shows percentage of CI that is in the null region (the ROPE)
#the null hypothesis is rejected or accepted if the percentage of the posterior 
#within the ROPE is smaller than to 2.5% or greater than 97.5%. 
#Desirable results are low proportions inside the ROPE (closer to zero the better).

outcome_tidy <- ebs_pop_final %>% 
  tidy_draws() %>% 
  rename(ate = `bs_stemperature_1`) #average temperature effect

# Find the proportion of posterior draws that are less than 0
outcome_tidy %>% 
  summarize(prop_lessthan_0 = sum(ate < 0) / n())
#99% chance that the average temperature effect is negative

#calculate the proportion of the posterior distribution that falls within ROPE
#set ROPE
rope <- 0.1 * sd(ebs.dat$Total_FA_Conc_WWT) # +/- 5.9

# Find the proportion of posterior draws that are less than 0
prop_outside <- outcome_tidy %>% 
  summarize(prop_outside_rope = 1 - sum(ate >= -5.9 & ate <= 5.9) / n())
#99.8% of the posterior distribution lies outside the ROPE
#pretty strong evidence of a large temperature effect 

#plot
ggplot(outcome_tidy, aes(x = ate)) +
  stat_halfeye(aes(fill_ramp = after_stat(x >= 5.9 | x <= -5.9)), 
               fill = "#4292C6", .width = c(0.95, 0.8)) +
  #scale_fill_ramp_manual(my_colors, guide = "none") +
  annotate(geom = "rect", xmin = -5.9, xmax = 5.9, ymin = -Inf, ymax = Inf, 
           fill = "red", alpha = 0.3) +
  labs(x = "Average temperature effect on energetic condition", y = NULL,
       caption = "Median shown with point; 80% and 95% credible intervals shown with black bars",
       y = "Posterior density") +
  theme_minimal()
#Given that we can't interpret credible intervals/estimates of smoothed 
  #coefficients, I don't think this is a valid approach 
#Also I think we'd need to standardize all response and predictor variables to obtain
  #comparable effect sizes 

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
plot(conditional_effects(mod6), ask = FALSE, points = TRUE, rug =TRUE)
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
loo2 <- loo(ebs_yrixn_final, save_psis = TRUE)
plot(loo2)

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
ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_yrixn_final))) #slight overdispersion
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_yrixn_final))

ppc_loo_pit_qq(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_yrixn_final),
               lw = weights(loo2$psis_object)) #hmmm..also not great

#Not going to tackle for this iteration, but it may be worthwhile to set more informative 
  #priors for this model. Diagnostic checks indicate that model isn't capturing extent of data 

################################
#Extract and plot conditional effects of yr*cpue and yr*temperature interaction

#temperature*year
conditions <- data.frame(year = c(2019, 2021, 2022, 2023))
ce1s_1 <- conditional_effects(ebs_yrixn_final, effects = "temperature",conditions = conditions, re_formula = NA,
                    probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(ebs_yrixn_final, effects = "temperature",conditions = conditions, re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_yrixn_final, effects = "temperature",conditions = conditions, re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$temperature
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temperature[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temperature[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temperature[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temperature[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(year), fill = ordered(year))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(year)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(year)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(year)), linewidth=1) +
  geom_rug(data = ebs.dat, aes(x = temperature, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") +  #raw data 
 theme_minimal() +
  labs(x = "Temperature", y = "Energetic Condition") +
  theme(legend.position="bottom") +
  theme(legend.position="none") +
   scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) -> temp

#density*year
conditions <- data.frame(year = c(2019, 2021, 2022, 2023))
ce1s_1 <- conditional_effects(ebs_yrixn_final, effects = "fourth.root.cpue",conditions = conditions, re_formula = NA,
                              probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(ebs_yrixn_final, effects = "fourth.root.cpue",conditions = conditions, re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(ebs_yrixn_final, effects = "fourth.root.cpue",conditions = conditions, re_formula = NA,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(year), fill = ordered(year))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(year)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(year)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(year)), linewidth = 1) +
  geom_rug(data = ebs.dat, aes(x = fourth.root.cpue, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") +  #raw data 
  theme_minimal() +
  labs(x = "Snow Crab Density", y = "Energetic Condition") +
  theme(legend.position="none") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) -> cpue

#Combine EBS plots 
(cpue + temp)  +  plot_annotation(tag_levels = 'a', title = "Eastern Bering Sea Snow Crab",
theme = theme(plot.title = element_text(hjust = 0.5))) -> ebs
                            
#Now run lines 740 - 809 in "NBS drivers models.R" to comine EBS and NBS plots
  #Using ggarrange instead of patchwork b/c both titles are retained!
ggarrange(ebs, nbs, nrow = 2, ncol = 1)

ggsave("./figures/yrinteractions_Fig6.png", plot=combine_plot,
       width = 6.5, height = 6.5, units = "in")

#faceted year x temperature plot
ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = .3, colour = NA) +
  geom_line(size=1) +
  geom_rug(data = ebs.dat, aes(x = temperature, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") +  #raw data 
  theme_minimal() +
  labs(x = "Temperature", y = "Energetic Condition") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~year)

#faceted year x density plot
ggplot(dat_ce, aes(x = effect1__, y = estimate__)) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90), alpha = .3, colour = NA) +
  geom_line(size=1) +
  geom_rug(data = nbs.dat, aes(x = fourth.root.cpue, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") +  #raw data 
  theme_minimal() +
  labs(x = "Snow Crab Density ", y = "Energetic Condition") +
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~year)

















