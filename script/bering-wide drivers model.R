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
library(DHARMa)
library(priorsense)
source("./script/stan_utils.R")

#load data
condition_master <- read.csv("./data/total_FA_master.csv")

#colors for plotting
my_colors <- c("#084594", "#9ECAE1", "#FCAE91", "#D55E00", "#4292C6")

my_colors3 <- c("#D55E00", "#084594")

#functions
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

#--------------------------------------------------------------------------------
#data wrangling
condition_master %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska"))) %>%  #add julian date 
  filter(!vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66"), 
         maturity != 1,
         sample_region > 0,
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

basin1_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) + 
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
#These models are very similar in out of sample predictive skill, but LOO favors 
  #model 2 with CPUE

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
#Temperature model has highest predictive capacity over base model and base +
  #CPUE model. Lets look at cpue + temp model with interaction

#------------------------------------------------------------------------------
#MODEL 4 BASE + TEMP x DENSITY x LME INTERACTION: default priors, truncated Gaussian, 
  #Note that while we could test a more simplistic cpue/temp additive model only before this, 
  #the interaction still allows use to test for fixed temp/cpue effects, and our a priori 
  #hypothesis based on Szuwalski et al is that we should expect an interaction, and that the 
  #strength of this interaction may differ between EBS and NBS

basin4_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      fourth.root.cpue*temperature*lme + (1 | region))   

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

#Interaction plot
conditions <- make_conditions(basin4, "lme")
conditional_effects(basin4, "fourth.root.cpue:temperature", conditions=conditions)

summary(basin4) 
bayes_R2(basin4) #R2 = 0.17
loo(basin4) 

# model comparison
loo(basin1, basin2, basin3, basin4, moment_match = TRUE)
#

#--------------------------------------------------------------------------------
#Full model Comparison (base model + base/cpue + base/temp + base/cpue*temp)

#LOO-CV
basin1 <- add_criterion(basin1, "loo")
basin2 <- add_criterion(basin2, "loo")
basin3 <- add_criterion(basin3, "loo")
basin4 <- add_criterion(basin4, "loo")
loo_compare(basin1, basin2, basin3, basin4, criterion = "loo") %>% print(simplify = F)
model.comp <- loo(basin1, basin2, basin3, basin4, moment_match = TRUE)

#and loo weights
model_weights(basin1, basin2, basin3, basin4, weights = "loo") %>% round(digits = 2)

#Table of Rsq Values 
rbind(bayes_R2(basin1), 
      bayes_R2(basin2),
      bayes_R2(basin3),
      bayes_R2(basin4)) %>%
  as_tibble() %>%
  mutate(model = c("basin1", "basin2", "basin3", "basin4"),
         r_square_posterior_mean = round(Estimate, digits = 2)) %>%
  select(model, r_square_posterior_mean) 

#weights 
loo1 <- loo(basin1)
loo2 <- loo(basin2)
loo3 <- loo(basin3)
loo4 <- loo(basin4)

loo_list <- list(loo1, loo2, loo3, loo4)

#Compute and compare Pseudo-BMA weights without Bayesian bootstrap, 
#Pseudo-BMA+ weights with Bayesian bootstrap, and Bayesian stacking weights
stacking_wts <- loo_model_weights(loo_list, method="stacking")
pbma_BB_wts <- loo_model_weights(loo_list, method = "pseudobma")
pbma_wts <- loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
round(cbind(stacking_wts, pbma_wts, pbma_BB_wts),2)

#Save model output 
tab_model(basin1, basin2, basin3, basin4)

forms <- data.frame(formula=c(as.character(basin1_formula)[1],
                              as.character(basin2_formula)[1],
                              as.character(basin3_formula)[1],
                              as.character(basin4_formula)[1]))

comp.out <- cbind(forms, model.comp$diffs[,1:2])

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
summary(basin_pop_final) #temp x cpue interaction differs between regions
bayes_R2(basin_pop_final) #r2 = .17
loo1 <- loo(basin_pop_final, save_psis = TRUE)
plot(loo1)

#Diagnostic Plots
plot(basin_pop_final, ask = FALSE)
plot(conditional_effects(basin_pop_final), ask = FALSE)
plot(conditional_smooths(basin_pop_final), ask = FALSE)
mcmc_plot(basin_pop_final, prob = 0.95)
mcmc_neff(neff_ratio(basin_pop_final)) #Effective sample size: All ratios > 0.1

#DHARMa residuals
  #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
  #https://frodriguezsanchez.net/post/using-dharma-to-check-bayesian-models-fitted-with-brms/#:~:text=Model%20checking%20with%20DHARMa,with%20other%20supported%20model%20types.
model.check <- createDHARMa(
  simulatedResponse = t(posterior_predict(basin_pop_final)),
  observedResponse = basin.dat$Total_FA_Conc_WWT,
  fittedPredictedResponse = apply(t(posterior_epred(basin_pop_final)), 1, mean),
  integerResponse = TRUE)

plot(model.check)
plot(model.check, form = basin.dat$region)
testDispersion(model.check)

#Interaction plot
conditions <- make_conditions(basin_pop_final, "lme")
conditional_effects(basin_pop_final, "fourth.root.cpue:temperature", conditions=conditions)

#Posterior Predictive Check Plots:
pp_check(basin_pop_final)
pp_check(basin_pop_final, type = "ecdf_overlay")
pp_check(basin_pop_final, type = "stat", stat = "mean")
pp_check(basin_pop_final, type = "stat", stat = "min")
pp_check(basin_pop_final, type = "stat", stat = "max")

#Marginal posterior predictive checks: Pit plots
ppc_pit_ecdf(pit=pit(y = basin.dat$Total_FA_Conc_WWT, yrep = posterior_predict(basin_pop_final))) #no overdisersion, looks good
ppc_intervals(y = basin.dat$Total_FA_Conc_WWT, yrep = posterior_predict(basin_pop_final))

ppc_loo_pit_qq(y = basin.dat$Total_FA_Conc_WWT, yrep = posterior_predict(basin_pop_final),
               lw = weights(loo1$psis_object)) #looks fairly uniform, no clear model mispecifications

#--------------------------------------------------------------------------------
#Extract and plot conditional effects of each predictor from best model (i.e. posterior distributions of conditional means)
  #conditioning on the mean for all other predictors, yr/site effects ignored

##NOTE: Final MS figures for Obj 2-3 are below following revisions

#Size effect plot for EBS/NBS
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

##Julian Day plot for EBS/NBS
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

#######
#Figure 4a: EBS temperature x density interaction

#Define temperature levels for interaction
temp <- list(temperature = c(0, 1, 2, 3))

#generate conditions for all levels of lme- use EBS level only for plot
conditions_all <- make_conditions(basin_pop_final, vars = "lme")
conditions_one_level <- conditions_all[conditions_all$lme == "EBS", ]

ce1s_1 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(temperature), fill = ordered(temperature))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(temperature)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(temperature)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(temperature)), linewidth = 1) +
  theme_minimal() +
  labs(x = "Snow Crab Density\n(Fourth root CPUE)", y = "Energetic Condition\n(mg FA/g WWT)",
       fill = expression("Temperature " ( degree~C)), color = expression("Temperature " ( degree~C))) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) + 
  ggtitle("Collapsing Eastern Bering Sea") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  theme(axis.title = element_text(size=10)) +
  theme(legend.title=element_text(size=9.5)) -> ebs_ixnplot

#####
#NBS temp x density interaction

#Define temperature levels for interaction
temp <- list(temperature = c(0, 1, 2, 3))

#generate conditions for all levels of lme- use NBS level only for plot
conditions_all <- make_conditions(basin_pop_final, vars = "lme")
conditions_one_level <- conditions_all[conditions_all$lme == "NBS", ]

ce1s_1 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:temperature",
                              int_conditions = temp, re_formula = NA,
                              conditions=conditions_one_level,
                              probs = c(0.1, 0.9))
dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(temperature), fill = ordered(temperature))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(temperature)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(temperature)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(temperature)), linewidth = 1) +
  theme_minimal() +
  labs(x = "Snow Crab Density\n(Fourth root CPUE)", y = "Energetic Condition\n(mg FA/g WWT)",
       fill = expression("Temperature " ( degree~C)), color = expression("Temperature " ( degree~C))) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) + 
  ggtitle("Non-collapsing Northern Bering Sea") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  theme(axis.title = element_text(size=10)) +
  theme(legend.title=element_text(size=9.5)) -> nbs_ixnplot

#####
#Temp x lme interaction 

ce1s_1 <- conditional_effects(basin_pop_final, effects = "temperature:lme",
                              re_formula = NA,
                              probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final, effects = "temperature:lme",
                              re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final, effects = "temperature:lme",
                              re_formula = NA,
                              probs = c(0.1, 0.9))

dat_ce <- ce1s_1$temperature
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$temperature[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$temperature[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$temperature[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$temperature[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(lme), fill = ordered(lme))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(lme)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(lme)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(lme)), linewidth = 1) +
  theme_minimal() +
  labs(x = expression("Temperature " ( degree~C)), y = "Energetic Condition\n(mg FA/g WWT)",
       fill = "", color = "") +
  scale_fill_manual(values = my_colors3) +
  scale_color_manual(values = my_colors3) + 
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  theme(axis.title = element_text(size=10)) +
  theme(legend.title=element_text(size=9.5)) -> temp_lme

#Density x lme interaction 

ce1s_1 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:lme",
                              re_formula = NA,
                              probs = c(0.025, 0.975))

## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:lme",
                              re_formula = NA,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final, effects = "fourth.root.cpue:lme",
                              re_formula = NA,
                              probs = c(0.1, 0.9))

dat_ce <- ce1s_1$fourth.root.cpue
dat_ce[["upper_95"]] <- dat_ce[["upper__"]]
dat_ce[["lower_95"]] <- dat_ce[["lower__"]]
dat_ce[["upper_90"]] <- ce1s_2$fourth.root.cpue[["upper__"]]
dat_ce[["lower_90"]] <- ce1s_2$fourth.root.cpue[["lower__"]]
dat_ce[["upper_80"]] <- ce1s_3$fourth.root.cpue[["upper__"]]
dat_ce[["lower_80"]] <- ce1s_3$fourth.root.cpue[["lower__"]]

ggplot(dat_ce, aes(x = effect1__, y = estimate__, color = ordered(lme), fill = ordered(lme))) +
  geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = ordered(lme)), alpha = .1, colour = NA) +
  geom_ribbon(aes(ymin = lower_90, ymax = upper_90, fill = ordered(lme)), alpha = .3, colour = NA) +
  geom_line(aes(color = ordered(lme)), linewidth = 1) +
  theme_minimal() +
  labs(x = "Snow Crab Density\n(Fourth root CPUE)", y = "Energetic Condition\n(mg FA/g WWT)",
       fill = "", color = "") +
  scale_fill_manual(values = my_colors3) +
  scale_color_manual(values = my_colors3) + 
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  theme(axis.title = element_text(size=10)) +
  theme(legend.title=element_text(size=9.5)) -> cpue_lme

#####
#Now just show NBS only temperature and density main effects by showing one level
  #of interaction only

#Temperature
#generate conditions for all levels of lme- use NBS level only for plot
conditions_all <- make_conditions(basin_pop_final, vars = "lme")
conditions_one_level <- conditions_all[conditions_all$lme == "NBS", ]

ce1s_1 <- conditional_effects(basin_pop_final, effects = "temperature",
                              re_formula = NA, conditions = conditions_one_level,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effects = "temperature", 
                              re_formula = NA, conditions = conditions_one_level,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effects = "temperature", 
                              re_formula = NA, conditions = conditions_one_level,
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
  theme(axis.title = element_text(size=10)) -> nbs_tempplot

##Snow Crab Density 
## 95% CI
ce1s_1 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", 
                              re_formula = NA, conditions = conditions_one_level,
                              probs = c(0.025, 0.975))
## 90% CI
ce1s_2 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", 
                              re_formula = NA, conditions = conditions_one_level,
                              probs = c(0.05, 0.95))
## 80% CI
ce1s_3 <- conditional_effects(basin_pop_final , effect = "fourth.root.cpue", 
                              re_formula = NA, conditions = conditions_one_level,
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
  theme(axis.title = element_text(size=10)) -> nbs_cpueplot

#---------------------------------------------------------------------------
#Combine and Save Plots to create Fig 4 and Supplementary Fig 

#Fig 4 temperature/density effects
ebsplot <- plot_spacer() + ebs_ixnplot + plot_spacer() + 
  plot_layout(widths = c(1,4.5,1))

nbsplot <- nbs_tempplot + nbs_cpueplot + 
  plot_annotation('Non-collapsing Northern Bering Sea',
                  theme=theme(plot.title=element_text(hjust=0.5, size=12)))

wrap_elements(ebsplot + plot_annotation(tag_levels='a')) / wrap_elements(nbsplot + plot_annotation(tag_levels = list(c('b','c'))))
ggsave("./figures/Fig4.png", height=7, width=6, unit="in")

#Supplementary Figure with size and DOY conditional effects 

(dayplot_basin + labs(y="Energetic Condition\n(mg FA/g WWT)"))  + sizeplot_basin +
  plot_annotation(tag_levels = 'a', 
                  theme = theme(plot.title = element_text(hjust = 0.5))) -> basin_supp
ggsave("./figures/Sup1.png", height=4, width=7, unit="in")
