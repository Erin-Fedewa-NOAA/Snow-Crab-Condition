#Based on reviewer comments: fitting a Bering-wide basin scale model (EBS + NBS combined) to 
  #investigate the effects of temperature and snow crab density on energetic condition 

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
#data wrangling
condition_master %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(!vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"), 
         maturity != 1,
         Total_FA_Conc_WWT > 0) %>%
  mutate(year = as.factor(year),
         lme = as.factor(lme),
         sex = as.factor(sex),
         region = as.factor(sample_region),
         station = as.factor(gis_station),
         temperature = as.numeric(gear_temperature),
         fourth.root.cpue = as.numeric(cpue^0.25),
         cpue = as.numeric(cpue),
         invert = as.numeric(total_benthic_cpue),
         julian = as.numeric(julian)) -> basin.dat

#Assess collinearity b/w covariates 
basin.dat %>%
  group_by(year, station) %>%
  summarise(temperature = mean(temperature),
            cpue = mean(cpue), 
            invert = mean(invert), 
            julian = mean(julian)) -> corr.dat

cor(corr.dat[,3:6]) 
corrplot(cor(corr.dat[,3:6]), method = 'number') 

#---------------------------------------------------------------------------------------
#MODEL 1 base model: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect (all nuisance sampling design covariates)

basin1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) + lme +
                      (1 | region)) 

basin1 <- brm(basin1_formula,
            data = basin.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(basin1, file = "./output/basin1.rds")
basin1 <- readRDS("./output/basin1.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin1$fit)
neff_lowest(basin1$fit)
rhat_highest(basin1$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(basin1, ask = FALSE)
plot(conditional_effects(basin1), ask = FALSE)
mcmc_plot(basin1, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(basin1)) #Effective sample size: All ratios > 0.1
pp_check(basin1) #not bad!

summary(basin1) #credible intervals for spline variance parameters (sds) don't include 0, let's keep smooths
bayes_R2(basin1) #R2 = 0.14
loo(basin1)

#----------------------------------------------------------------------------------
#MODEL 2 BASE + CRAB CPUE BY LME: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + crab cpue x lme interaction 

basin2_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(fourth.root.cpue, k = 3, by = lme) + (1 | region))  

basin2 <- brm(basin2_formula,
            data = basin.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))


#Save output
saveRDS(basin2, file = "./output/basin2.rds")
basin2 <- readRDS("./output/basin2.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin2$fit)
neff_lowest(basin2$fit)
rhat_highest(basin2$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(basin2, ask = FALSE)
plot(conditional_effects(basin2), ask = FALSE)
mcmc_plot(basin2, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(basin2)) #Effective sample size: All ratios > 0.1
pp_check(basin2)

summary(basin2) 
bayes_R2(basin2) #R2 = 0.14
loo(basin2) 

# model comparison
loo(basin1, basin2, moment_match = TRUE) 
#Really no difference in these two models. We'll keep cpue in as a way to test hypotheses
#in full model 

#--------------------------------------------------------------------------------
#MODEL 3 BASE + TEMP: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature fixed effect 

basin3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(temperature, k=3, by=lme) + (1 | region))  

basin3 <- brm(basin3_formula,
            data = basin.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(basin3, file = "./output/basin3.rds")
basin3 <- readRDS("./output/basin3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin3$fit)
neff_lowest(basin3$fit)
rhat_highest(basin3$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(basin3, ask = FALSE)
plot(conditional_effects(basin3), ask = FALSE)
mcmc_plot(basin3, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(basin3)) #Effective sample size: All ratios > 0.1
pp_check(basin3)

summary(basin3) 
bayes_R2(basin3) #R2 = 0.17
loo(basin3) -> d
plot(d)

# model comparison
loo(basin1, basin2, basin3, moment_match = TRUE)
#Temperature model has highest predictive capacity 

#------------------------------------------------------------------------------
#MODEL 4 BASE + TEMP + DENSITY: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature and density fixed effect 

basin4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      s(temperature, k=3, by=lme) + s(fourth.root.cpue, k=3, by=lme) +
                      (1 | region))  

basin4 <- brm(basin4_formula,
            data = basin.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(basin4, file = "./output/basin4.rds")
basin4 <- readRDS("./output/basin4.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin4$fit)
neff_lowest(basin4$fit)
rhat_highest(basin4$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(basin4, ask = FALSE)
plot(conditional_effects(basin4), ask = FALSE)
mcmc_plot(basin4, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(basin4)) #Effective sample size: All ratios > 0.1
pp_check(basin4)

summary(basin4) 
bayes_R2(basin4) #R2 = 0.18
loo(basin4) 

# model comparison
loo(basin1, basin2, basin3, basin4, moment_match = TRUE)
#So seems that full additive basin4 has highest predictive capacity 
#and is a substantial improvement over basin1 and basin2 but addition of cpue 
#doesnt improve much 

#------------------------------------------------------------------------------
#MODEL 5 BASE + TEMP x DENSITY x LME INTERACTION: default priors, truncated Gaussian, 
#Covariates: crab size/julian day/random effect + temperature x density x lme linear interaction

basin5_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      fourth.root.cpue*temperature*lme + (1 | region))   

basin5 <- brm(basin5_formula,
            data = basin.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500, warmup = 1000,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(basin5, file = "./output/basin5.rds")
basin5 <- readRDS("./output/basin5.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin5$fit)
neff_lowest(basin5$fit)
rhat_highest(basin5$fit) #Potential scale reduction: All rhats < 1.1

#Diagnostic Plots
plot(basin5, ask = FALSE)
plot(conditional_effects(basin5), ask = FALSE)
mcmc_plot(basin5, type = "areas", prob = 0.95)
mcmc_neff(neff_ratio(basin5)) #Effective sample size: All ratios > 0.1
pp_check(basin5)

#Interaction plot
conditions <- make_conditions(basin5, "lme")
conditional_effects(basin5, "fourth.root.cpue:temperature", conditions=conditions)

summary(basin5) 
bayes_R2(basin5) #R2 = 0.17
loo(basin5) 

# model comparison
loo(basin1, basin2, basin3, basin4, basin5, moment_match = TRUE)
#Looks like interaction did not improve predictive capacity from main effects model

#--------------------------------------------------------------------------------
#Full model Comparison (base model + base/cpue + base/cpue/temp + base/cpue*temp)

#LOO-CV
basin1 <- add_criterion(basin1, "loo")
basin3 <- add_criterion(basin3, "loo")
basin4 <- add_criterion(basin4, "loo")
basin5 <- add_criterion(basin5, "loo")
loo_compare(basin1, basin3, basin4, basin5, criterion = "loo") %>% print(simplify = F)
model.comp <- loo(basin1, basin3, basin4, basin5, moment_match = TRUE)

#and loo weights
model_weights(basin1, basin3, basin4, basin5, weights = "loo") %>% round(digits = 2)

#Table of Rsq Values 
rbind(bayes_R2(basin1), 
      bayes_R2(basin3),
      bayes_R2(basin4),
      bayes_R2(basin5)) %>%
  as_tibble() %>%
  mutate(model = c("basin1", "basin3", "basin4", "basin5"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#weights 
loo1 <- loo(basin1)
loo3 <- loo(basin3)
loo4 <- loo(basin4)
loo5 <- loo(basin5)

loo_list <- list(loo1, loo3, loo4, loo5)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
#Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)

#Save model output 
tab_model(basin1, basin2, basin3, basin4, basin5)

forms <- data.frame(formula=c(as.character(basin1_formula)[1],
                              as.character(basin2_formula)[1],
                              as.character(basin3_formula)[1],
                              as.character(basin4_formula)[1],
                              as.character(basin5_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])
write.csv(comp.out, "./output/basin_pop_model_comp.csv")

#--------------------------------------------------------------------------
#FINAL MODEL:  Run basin4 model with 10,000 iterations and set seed for reproducibility 
basin_pop_final <- brm(basin4_formula,
                     data = basin.dat,
                     family = gaussian,
                     cores = 4, chains = 4, iter = 10000, warmup = 1000,
                     save_pars = save_pars(all = TRUE), seed = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(basin_pop_final, file = "./output/basin_pop_final.rds")
basin_pop_final <- readRDS("./output/basin_pop_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(basin_pop_final$fit)
neff_lowest(basin_pop_final$fit)
rhat_highest(basin_pop_final$fit)
summary(basin_pop_final)
bayes_R2(basin_pop_final) #r2 = .11
loo1 <- loo(basin_pop_final, save_psis = TRUE)
plot(loo1)

#Diagnostic Plots
plot(basin_pop_final, ask = FALSE)
plot(conditional_effects(basin_pop_final), ask = FALSE)
mcmc_plot(basin_pop_final, prob = 0.95)
mcmc_neff(neff_ratio(basin_pop_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(basin_pop_final)
pp_check(basin_pop_final, type = "ecdf_overlay")
pp_check(basin_pop_final, type = "stat", stat = "mean")
pp_check(basin_pop_final, type = "stat", stat = "min")
pp_check(basin_pop_final, type = "stat", stat = "max")

#Marginal posterior predictive check: Pit plot
ppc_loo_pit_qq(y = basin.dat$Total_FA_Conc_WWT, yrep = posterior_predict(basin_pop_final),
               lw = weights(loo1$psis_object)) 

#--------------------------------------------------------------------------------
#Extract and plot conditional effects of each predictor from best model (i.e. posterior distributions of conditional means)
#conditioning on the mean for all other predictors, yr/site effects ignored 

#Size effect plot 
#Need to save settings from conditional effects as an object to plot in ggplot
## 95% CI
ce1s_1 <- conditional_effects(basin_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effect = "cw", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effect = "cw", re_formula = NA,
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
  geom_rug(data = basin.dat, aes(x = cw, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data 
  labs(x = "Carapace Width (mm)", y = "") +
  theme_minimal() +
  ylim(0,215) +
  theme(axis.text=element_text(size=8)) +
  theme(axis.title = element_text(size=10)) -> sizeplot_basin

##Julian Day
## 95% CI
ce1s_1 <- conditional_effects(basin_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effect = "julian", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effect = "julian", re_formula = NA,
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
  geom_rug(data = basin.dat, aes(x = julian, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Day of Year", y = "") +
  theme_minimal() +
  ylim(0,215) +
  theme(axis.title = element_text(size=10)) -> dayplot_basin

##Snow Crab Density 
## 95% CI
ce1s_1 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", re_formula = NA,
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
  geom_rug(data = basin.dat, aes(x = fourth.root.cpue, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = "Snow Crab Density\n(Fourth root CPUE)", y = "") +
  theme_minimal() +
  ylim(0,215) +
  theme(axis.title = element_text(size=10)) -> cpueplot_basin

##Temperature 
## 95% CI
ce1s_1 <- conditional_effects(basin_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effect = "temperature", re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effect = "temperature", re_formula = NA,
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
  geom_rug(data = basin.dat, aes(x = temperature, y = Total_FA_Conc_WWT), 
           colour = "grey80", linewidth = .5, sides="b", alpha=.7, position = "jitter") + #raw data
  labs(x = expression("Temperature " ( degree~C)), y = "Energetic Condition\n(mg FA/g WWT)") +
  theme_minimal() +
  ylim(0, 215) +
  theme(axis.title = element_text(size=10)) -> tempplot_basin

#See EBS drivers script for figure compilation 
