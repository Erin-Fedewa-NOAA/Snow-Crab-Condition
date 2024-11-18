#Estimate condition for each year/region (EBS/NBS) after controlling for the effect of size due 
  #to difference in ontogeny/sex, seasonality due to sampling design of survey, and spatial 
  #variation within the EBS/NBS (region group-level effect)

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
#2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
#used as response variable in all further analyses. 

# Author: EJF

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
library(hrbrthemes)
library(bayestestR)
source("./script/stan_utils.R")

#load data
condition_master <- read.csv("./data/total_FA_master.csv")

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

#colors
my_colors <- c("#D55E00","#9ECAE1", "#4292C6", "#084594")
my_colors2 <- RColorBrewer::brewer.pal(7, "GnBu")[c(5,6)]


#############################################
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
         cpue = as.numeric(cpue),
         invert = as.numeric(total_benthic_cpue),
         julian = as.numeric(julian)) -> ebs.dat 

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
         cpue = as.numeric(cpue),
         invert = as.numeric(total_benthic_cpue),
         julian = as.numeric(julian)) -> nbs.dat

################################################
#EBS ANNUAL MEANS

ebs_annual_final_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                      year + (1 | region))  

ebs_annual_final <- brm(ebs_annual_final_formula,
                     data = ebs.dat,
                     family = gaussian,
                     cores = 4, chains = 4, iter = 10000, warmup = 1000,
                     save_pars = save_pars(all = TRUE), seed = 3,
                     control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(ebs_annual_final, file = "./output/ebs_annual_final.rds")
ebs_annual_final <- readRDS("./output/ebs_annual_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(ebs_annual_final$fit)
neff_lowest(ebs_annual_final$fit)
rhat_highest(ebs_annual_final$fit)
summary(ebs_annual_final) #dramatically lower condition in 2019
bayes_R2(ebs_annual_final) #r2 = .38 - including year in the base model greatly improves fit
loo(ebs_annual_final, moment_match = T)

#Diagnostic Plots
plot(ebs_annual_final, ask = FALSE)
plot(conditional_effects(ebs_annual_final), ask = FALSE)
mcmc_plot(ebs_annual_final, prob = 0.95)
mcmc_neff(neff_ratio(ebs_annual_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(ebs_annual_final)
pp_check(ebs_annual_final, type = "ecdf_overlay")
pp_check(ebs_annual_final, type = "stat", stat = "mean")
pp_check(ebs_annual_final, type = "stat", stat = "min")
pp_check(ebs_annual_final, type = "stat", stat = "max")

#Pit plots
ppc_pit_ecdf(pit=pit(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_annual_final))) #no overdisersion, looks good
ppc_intervals(y = ebs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(ebs_annual_final))

################################
#Extract conditional effect of year for plot

conditional_effects(ebs_annual_final, effect = "year")

ce1s_1 <- conditional_effects(ebs_annual_final, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(lme = as.factor("Eastern Bering Sea"))-> year_ebs

#And some playing around with marginal effects/ROPE

#Average marginal effect of year (difference across yrs while holding cw and julian day constant)
years_ame <- ebs_annual_final %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NULL) %>% 
  gather_emmeans_draws()

years_ame %>%
  median_hdi()

#ebs marginal effects plot 
ggplot(years_ame,aes(x = .value, fill=ordered(year))) +
  stat_halfeye(slab_alpha = 0.7) +
  labs(x = "Average marginal effect",
       y = "Density") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  ggtitle("Eastern Bering Sea") +
  theme(plot.title = element_text(hjust = 0.5)) -> marg_ebs

#And a hypothetical example of how you might annotate a domain-specific ROPE
  #this would be a neat way to illustrate lab-derived FA thresholds! 
ggplot(years_ame, aes(y = year, x = .value, fill = after_stat(abs(x) < 50))) +
  stat_halfeye() +
  geom_vline(xintercept = c(0, 50), linetype = "dashed") +
  scale_fill_manual(values = c("gray80", "skyblue")) +
  theme_minimal() + 
  labs(x = "Average Marginal Effect on Energetic Condition", y="")

##########################
#Run NBS model

nbs_annual_final_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 3) + s(julian, k = 3) +
                                  year + (1 | region))  

nbs_annual_final <- brm(nbs_annual_final_formula,
                        data = nbs.dat,
                        family = gaussian,
                        cores = 4, chains = 4, iter = 10000, warmup = 1000,
                        save_pars = save_pars(all = TRUE), seed = 3,
                        control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save model output 
saveRDS(nbs_annual_final, file = "./output/nbs_annual_final.rds")
nbs_annual_final <- readRDS("./output/nbs_annual_final.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(nbs_annual_final$fit)
neff_lowest(nbs_annual_final$fit)
rhat_highest(nbs_annual_final$fit)
summary(nbs_annual_final)
bayes_R2(nbs_annual_final) #r2 = .11 
loo(nbs_annual_final)

#Diagnostic Plots
plot(nbs_annual_final, ask = FALSE)
plot(conditional_effects(nbs_annual_final), ask = FALSE)
mcmc_plot(nbs_annual_final, prob = 0.95)
mcmc_neff(neff_ratio(nbs_annual_final)) #Effective sample size: All ratios > 0.1

#Posterior Predictive Check Plots:
pp_check(nbs_annual_final) #looks much better than ebs dataset!
pp_check(nbs_annual_final, type = "ecdf_overlay")
pp_check(nbs_annual_final, type = "stat", stat = "mean")
pp_check(nbs_annual_final, type = "stat", stat = "min")
pp_check(nbs_annual_final, type = "stat", stat = "max")

#Pit plots
ppc_pit_ecdf(pit=pit(y = nbs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(nbs_annual_final))) #no overdisersion, looks good
ppc_intervals(y = nbs.dat$Total_FA_Conc_WWT, yrep = posterior_predict(nbs_annual_final)) 

######################################
#Conditional Effect for year 
conditional_effects(nbs_annual_final, effect = "year")

ce1s_1_nbs <- conditional_effects(nbs_annual_final, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1_nbs$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(lme = as.factor("Northern Bering Sea")) -> year_nbs

#Combine EBS and NBS plots for Fig 3 in ms 
year_ebs %>%
  full_join(year_nbs) -> dat2
#now plot
ggplot(dat2, aes(year, estimate__,)) +
  geom_bar(aes(fill = ordered(year)), stat='identity', size=3) +
  geom_errorbar(aes(year, ymin=lower__, ymax=upper__), width=0.3, size=0.5, color = "grey60") +
  ylab("Energetic Condition (Total FA per DWT)") + xlab("") +
  
  scale_fill_manual(values=my_colors) +
  facet_wrap(~lme) +
  geom_vline(data = subset(dat2, lme == "Eastern Bering Sea"), aes(xintercept = 1.5), linetype="dashed") +
  geom_text(data = subset(dat2, lme == "Eastern Bering Sea"), aes(x = 1, y=185, label = "Mid-collapse"),
            size = 2.5, color = "#D55E00") +
  geom_text(data = subset(dat2, lme == "Eastern Bering Sea"), aes(x = 3, y=185, label = "Post-collapse"),
            size = 2.5, color = "#084594") +
  theme_ipsum(axis_title_just = "cc", axis_title_size = 13, axis_text_size =11) +
  theme(legend.position="none") +
  theme(panel.grid.major.x = element_blank()) + 
  theme(axis.title.y = element_text(size = 12)) 

ggsave("./figures/Fig3.png", dpi=300, width = 7, height = 5, units = "in")

#Average marginal effect of year (difference across yrs while holding cw and julian day constant)
years_ame_nbs <- nbs_annual_final %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NULL) %>% 
  gather_emmeans_draws()

#nbs marginal effects plot 
ggplot(years_ame_nbs,aes(x = .value, fill=ordered(year))) +
  stat_halfeye(slab_alpha = 0.7) +
  labs(x = "Average marginal effect",
       y = "Density") +
  scale_fill_manual(values = my_colors) +
  theme_minimal() + 
  theme(legend.position="bottom") +
  theme(legend.title=element_blank()) + 
  ggtitle("Northern Bering Sea") +
  theme(plot.title = element_text(hjust = 0.5)) -> marg_nbs

#combine plots
marg_ebs + marg_nbs + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("./figures/FigSupp.png", dpi=300, width = 6.5, height = 4.5, units = "in")

