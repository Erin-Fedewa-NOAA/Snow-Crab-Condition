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
  geom_histogram(aes(log(cpue + 1)))

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
  ggplot(aes(log(Total_FA_Conc_WWT))) + 
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
pp_check(model_normal_trunc, ndraws = 1e2) + pp_check(model_normal_log, ndraws = 1e2) + pp_check(model_skew_trunc, ndraws = 1e2) +
  pp_check(model_skew_log, ndraws = 1e2)+ pp_check(model_gamma, ndraws = 1e2) + pp_check(model_log, ndraws = 1e2) 
#truncated skew normal captures mean and variance best, none of the models picking up apparent bimodality of data 

# posterior predictive checking - boxplots
pp_check(model_normal_trunc, type = "boxplot", ndraws = 20) + pp_check(model_normal_log, type = "boxplot", ndraws = 20) + pp_check(model_skew_trunc, type = "boxplot", ndraws = 20) +
  pp_check(model_skew_log, type = "boxplot", ndraws = 20) + pp_check(model_gamma, type = "boxplot", ndraws = 20) + 
  pp_check(model_log, type = "boxplot", ndraws = 20) 

#let's look at the distribution of minimum values for posterior distributions vrs data
pp_check(model_normal_trunc, type = "stat", stat = "min") + pp_check(model_normal_log, type = "stat", stat = "min") + pp_check(model_skew_trunc, type = "stat", stat = "min") +
  pp_check(model_skew_log, type = "stat", stat = "min") + pp_check(model_gamma, type = "stat", stat = "min") + 
  pp_check(model_log, type = "stat", stat = "min")

#now maximum values
pp_check(model_normal_trunc, type = "stat", stat = "max") + pp_check(model_normal_log, type = "stat", stat = "max") + pp_check(model_skew_trunc, type = "stat", stat = "max") +
  pp_check(model_skew_log, type = "stat", stat = "max") + pp_check(model_gamma, type = "stat", stat = "max") + 
  pp_check(model_log, type = "stat", stat = "max")
#lognormal is way overshooting max 

#and means
pp_check(model_normal_trunc, type = "stat", stat = "mean") + pp_check(model_normal_log, type = "stat", stat = "mean") + pp_check(model_skew_trunc, type = "stat", stat = "mean") +
  pp_check(model_skew_log, type = "stat", stat = "mean") + pp_check(model_gamma, type = "stat", stat = "mean") + 
  pp_check(model_log, type = "stat", stat = "mean")
  
#Which models are providing best predictions?
model_normal_trunc <- add_criterion(model_normal_trunc, "waic")
model_skew_trunc <- add_criterion(model_skew_trunc, "waic")
model_gamma <- add_criterion(model_gamma, "waic")
model_log <- add_criterion(model_log, "waic")
  loo_compare(model_normal_trunc, model_skew_trunc, model_gamma, model_log, criterion = "waic")
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
#Model runs not shown here, but group-level effects structure was explored. Due to the high number of stations\
  #containing only 1 crab, 1|station and nested 1|region/station models had convergence issues. 
#We'll go with 1|region to at least attempt to account for the repeat sampling design

#####
#MODEL 1 BASE MODEL: default priors, truncated Gaussian, 
    #Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)

mod1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      (1 | region))  

mod1 <- brm(mod1_formula,
            data = ebs.dat,
            family = gaussian,
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
bayes_R2(mod1) #R2 = 0.2
loo(mod1) -> a
plot(a)

#####
#MODEL 2 BASE MODEL + INVERT: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year interaction

mod2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + (1 | region))  

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
bayes_R2(mod2) #R2 = 0.29
loo(mod2) -> b
plot(b)

# model comparison
loo(mod1, mod2, moment_match = TRUE) 
#Looks like benthic invert density increases predictive skill 

#####
#MODEL 3 BASE MODEL + INVERT + CRAB: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year + crab*year interaction

mod3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + s(fourth.root.cpue, k=4, by = year) +
                      (1 | region))  

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
plot(conditional_smooths(mod3), ask = FALSE)
mcmc_plot(mod3, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(mod3)) #Effective sample size: All ratios > 0.1
pp_check(mod3)

summary(mod3) 
bayes_R2(mod3) #R2 = 0.45
loo(mod3) -> c
plot(c)

# model comparison
loo(mod1, mod2, mod3, moment_match = TRUE) 
#Snow crab CPUE also increases predictive capacity

#####
#MODEL 4 BASE MODEL + INVERT + CRAB + TEMP: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + invert*year + crab*year + temp*year interaction 

mod4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + s(julian, k = 4) +
                      s(fourth.root.invert, k=4, by = year) + s(fourth.root.cpue, k=4, by = year) +
                      s(temperature, k=4, by = year) + (1 | region))  

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
bayes_R2(mod4) #R2 = 0.49
loo(mod4) -> d
plot(d)

# model comparison
loo(mod1, mod2, mod3, mod4, moment_match = TRUE)
#So seems that full model has highest predictive capacity 
  
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



