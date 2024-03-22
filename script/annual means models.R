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

loo(mod3, moment_match = T)
loo(mod3.5, moment_match = T)
loo_compare(mod3, mod3.5)

m3 <- add_criterion(mod3, "waic")
m3.5 <- add_criterion(mod3.5, "waic")
loo_compare(mod3, mod3.5, criterion = "waic")

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