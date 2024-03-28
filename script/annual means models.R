#Predict condition for each year/region (EBS/NBS)

##NOTE: WWT:DWT ratios appear to be affected by difference in sampling methods in 
#2019. B/c total FA per WWT were not subject to the WWT:DWT discrepancy, it will be 
#used as response variable in all further analyses. 

mod3_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ s(cw, k = 4) + year + s(julian, k = 4) + (1 | region/station)) 

#seperate models for EBS and NBS, plot condition effect of year (so accounting for crab
  #size, julian date and )

#controlling for the effect of size due to difference in ontogeny/sex, seasonality 
  #due to sampling design of survey, and spatial variation within the EBS

mod3 <- brm(mod3_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(mod3, file = "./output/mod3.rds")
mod3 <- readRDS("./output/mod3.rds")

#MCMC convergence diagnostics 
check_hmc_diagnostics(mod2$fit)
neff_lowest(mod1$fit)
rhat_highest(mod1$fit)
pp_check(mod3, ndraws = 1000) + labs(title = str_glue("Posterior predictive checks for mod2"))
summary(mod2)
tidy(mod2)
bayes_R2(mod3)

plot(marginal_effects(mod3),points=T) #does it fit the data?
marginal_effects(mod3)
conditional_effects(mod3)

#better model if we drop crab size?
mod3.5_formula <-  bf(Total_FA_Conc_WWT | trunc(lb = 0) ~ year + s(julian, k = 4) + (1 | region/station)) 

mod3.5 <- brm(mod3.5_formula,
            data = ebs.dat,
            family = gaussian,
            cores = 4, chains = 4, iter = 2500,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

saveRDS(mod3.5, file = "./output/mod3.5.rds")
mod3.5 <- readRDS("./output/mod3.5.rds")

pp_check(mod3.5, ndraws = 1000) + labs(title = str_glue("Posterior predictive checks for mod2"))
summary(mod3.5)
bayes_R2(mod3.5)

plot(marginal_effects(mod3.5),points=T) #does it fit the data?
marginal_effects(mod3)
conditional_effects(mod3.5)

a <- loo(mod3) 
loo(mod3.5)
loo_compare(mod3, mod3.5)
k3 <- kfold(mod3, K=10)
k3.5 <- kfold(mod3.5, K=10)
loo_compare(mod3, mod3.5)

plot(a)

#Hmmm I think we're back to square one- region/station model is overfit 

##########################
#Run NBS model
#now extract annual conditional effect and plot, merge plots 


#Conditional Effect for year 
conditional_effects(hurdle1_tanner, effect = "year")

ce1s_1 <- conditional_effects(hurdle1_tanner, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Tanner crab") -> year_tanner

#######################
#spatial effects - too small of sample sizes to look by region- prob not of interest 

#Conditional Effect for year 
conditional_effects(hurdle1_snow, effect = "year")

ce1s_1 <- conditional_effects(hurdle1_snow, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Snow crab") -> year_snow

#Average marginal effect of year 
years_ame <- hurdle1_snow %>% 
  emmeans(~ year,
          var = "year",
          epred = TRUE, re_formula = NA) %>% 
  gather_emmeans_draws()

years_ame %>%
  median_hdi()

ggplot(years_ame,aes(x = .value, fill=year)) +
  stat_halfeye(slab_alpha = 0.75) +
  labs(x = "Average marginal effect",
       y = "Density") +
  theme_bw()

#Marginal effects: effect of year in the hurdling process 

#Combine hu_year term(s) with hurdle intercept and transform
hurdle_intercept <- tidy(hurdle1_snow) %>% 
  filter(term == "hu_(Intercept)") %>%  
  pull(estimate)

hurdle_lifeexp <- tidy(hurdle1_snow) %>%  
  filter(term == "hu_year2017") %>%  
  pull(estimate)

plogis(hurdle_intercept + hurdle_lifeexp) - plogis(hurdle_intercept)
#The probability of seeing 0% prevalence in 2017 decreased by 38% from 2015

#Hurdle Model 1: Ignoring spatial variability, year effect only (Obj 2 analysis)
hurdle1_formula <- bf(
  #mu, mean part of formula
  Prevalance ~ year,
  #alpha, zero inflation part
  hu ~ year) 

hurdle1_tanner <- brm(hurdle1_formula,
                      data = prev.tanner,
                      family = hurdle_lognormal(),
                      cores = 4, chains = 4, iter = 2500,
                      save_pars = save_pars(all = TRUE),
                      control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(hurdle1_tanner, file = "./output/hurdle1_tanner.rds")
hurdle1_tanner <- readRDS("./output/hurdle1_tanner.rds")

tidy(hurdle1_tanner)
pp_check(hurdle1_tanner) 
#both zero and non zero processes incorporated into the posterior distribution

#MCMC convergence diagnostics 
check_hmc_diagnostics(hurdle1_tanner$fit)
neff_lowest(hurdle1_tanner$fit)
rhat_highest(hurdle1_tanner$fit)
summary(hurdle1_tanner)
bayes_R2(hurdle1_tanner)

#Diagnostic Plots
plot(hurdle1_tanner, ask = FALSE)
conditional_effects(hurdle1_tanner)
mcmc_plot(hurdle1_tanner, prob = 0.95)
mcmc_plot(hurdle1_tanner, transformations = "inv_logit_scaled")

#Conditional Effect for year 
conditional_effects(hurdle1_tanner, effect = "year")

ce1s_1 <- conditional_effects(hurdle1_tanner, effect = "year", re_formula = NA,
                              probs = c(0.025, 0.975)) 
ce1s_1$year %>%
  dplyr::select(year, estimate__, lower__, upper__) %>%
  mutate(species = "Tanner crab") -> year_tanner

######################################
#Combine snow and tanner plots for Fig 5 in ms 
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
  ylab("Prevalance (%)") + xlab("") +
  scale_colour_manual(values = new_colors) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(legend.title= element_blank())
ggsave("./figs/annual_hurdle.png", dpi=300)