#Objective 1: Assess the relationship b/w total fatty acid concentration and hepato WWT:DWT
#Exploratory script assessing morphometric condition metrics relative to total FA/hepato WWT:DWT

#NOTE: We're removing 2019 from Obj 1 validation models b/c methods differed (i.e. hepatos were dissected
  #in the lab after freezing whole crab, and likely affects WWT:DWT ratio due to water loss)

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
library(mgcv)
library(viridis)
source("./script/stan_utils.R")

#load data 
condition_master <- read.csv("./data/total_FA_master.csv")

#plotting
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[c(3,5,6,7)]
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
options(scipen=999)

####################################
#Hepato %DWT condition metric 

#data wrangling 
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),#likely tanners
         maturity != 1,
         year > 2019) -> new.dat

#Let's look at distributions
condition_master %>%
  ggplot() +
  geom_histogram(aes(Perc_DWT))

condition_master %>%
  ggplot() +
  geom_histogram(aes(Total_FA_Conc_DWT))

#% Plot: DWT vrs total FA - faceted 
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA_Conc_DWT)) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA (mg/g DWT)") +
  facet_wrap(~year)

#% Plot: DWT vrs total FA concentration 
new.dat %>%
  ggplot(aes(Perc_DWT, Total_FA_Conc_DWT)) +
  geom_point(aes(color=factor(year))) +
  theme_bw() + 
  geom_smooth(method = "lm", colour="black", level = 0.95) +
  labs(x= "% DWT in Hepatopancreas", y = "Total FA per DWT (mg FA/g WWT)") +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values=cbPalette)

############################
#Bayesian regression model (FA ~ %DWT)
condition_1 <- brm(data = new.dat,
            family = gaussian(), #Student's t distribution- more robust to outliers
            Total_FA_Conc_DWT ~ 1 + Perc_DWT,
            seed=1,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(condition_1, file = "./output/condition_1.rds")
condition_1 <- readRDS("./output/condition_1.rds")

summary(condition_1)
bayes_R2(condition_1) #.64
posterior_summary(condition_1)
pairs(condition_1)
loo_c1 <- loo(condition_1)
pareto_k_table(loo_c1)

#Diagnostic Plots
plot(condition_1, ask = FALSE)

#Predictions
dwt_seq <- tibble(Perc_DWT = seq(from = 0, to = 60, by =1))

mu <- fitted(condition_1, newdata = dwt_seq) %>%
  as_tibble() %>%
  bind_cols(dwt_seq)

#Plot model fit - Fig 5 for Manuscript 
new.dat %>%
  ggplot(aes(x = Perc_DWT, y = Total_FA_Conc_DWT, color=year)) +
  #geom_abline(intercept = fixef(condition_1)[1], 
             # slope     = fixef(condition_1)[2],
             # size = .8, color = "black") +
  geom_smooth(data = mu, aes(y=Estimate, ymin= Q2.5, ymax= Q97.5, color=year), 
              stat = "identity", fill = "grey70", color = "black", alpha = 0.5) +
  geom_point(aes(color=as.factor(year)), size = 1) +
  theme(panel.grid = element_blank()) +
  theme_minimal() + 
  labs(x= "% DWT in Hepatopancreas", y = "Total FA per DWT (mg FA/g DWT)") +
  theme(legend.title=element_blank()) +
  scale_colour_manual(values=my_colors) +
  coord_cartesian(xlim = range(new.dat$Perc_DWT),
                  ylim = c(0,775)) +
  #make legend points bigger
  guides(colour = guide_legend(override.aes = list(size=4)))
ggsave("./figures/Fig5.png", height = 4, width = 5, units = "in", dpi = 300)

############################
#Exploring Morphometric condition metrics

#First, let's look at weight at size in survey data 
ebs_haul <- read.csv("./data/crabhaul_opilio.csv")
nbs_haul <- read.csv("./data/crabhaul_opilio_nbs.csv")

#Wrangle ebs data -   
ebs_haul %>%
  rename_with(tolower) %>%
  mutate(julian=yday(parse_date_time(start_date, "mdy", "US/Alaska")),  #add julian date
          year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  filter(year %in% c(2011, 2015, 2017:2024),
         haul_type==3, 
         weight != "NA", 
        weight != 0, #weird- one 0 observations...
         width != "NA") %>%
  mutate(log_weight = as.numeric(log(weight)),
         log_cw = as.numeric(log(width)), 
         year = as.factor(year),
         weight = as.numeric(weight)) -> dat

#sample sizes of all weight data
dat %>%
  group_by(year) %>%
  count()

#plot full dataset
dat %>%
  ggplot(aes(log_cw, log_weight,color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) 
#visually fits appear very similar across years 

#plot low sample size years only
dat %>%
  filter(year %in% c(2018, 2021, 2023)) %>%
  ggplot(aes(log_cw, log_weight, color=year)) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) 
#snow crab weights were not systematically collected in these years, and size range 
  #coverage of opportunistic samples isn't great for 2021 and 2023- let's drop  

dat %>% filter(year %in% c(2011, 2015, 2017:2019)) -> precollapse_dat
dat %>% filter(year %in% c(2011, 2015, 2017:2019, 2022, 2024)) -> full_dat

#plot by bottom temp-similar to Fig 49 in Science paper appendix
full_dat %>%
  ggplot(aes(width, weight, color=gear_temperature)) +
  geom_point() +
  theme_bw() + 
  scale_color_distiller(palette = "Spectral") +
  geom_smooth(method = "gam", se = FALSE, color="black", linewidth=.5) +
  facet_wrap(~year)
#seems the temperature covariate primarily accounts for weight at size 
#differences due to size/stage specific thermal preferences (i.e. positive effect of 
#gear temperature on weight simply b/c larger crab at warmer temps), but we'll keep in

#2017 vrs 2018 only
dat %>%
  filter(year %in% c(2017,2018)) %>%
  ggplot(aes(width, weight, color=gear_temperature, group=year)) +
  geom_point() +
  theme_bw() + 
  scale_color_distiller(palette = "Spectral") +
  geom_smooth(method = "gam", se = FALSE, color="black", linewidth=.8) 

#fit gam using Szuwalski et al approach 
mod1 <- gam(weight ~ s(width) + year + s(gear_temperature, k=3), data = precollapse_dat)
summary(mod1) #negative coefficient for 2018, but very small effect - not sure 
#how data was filtered down to n=27 in 2018?
plot.gam(mod1, se=TRUE, shade=TRUE, all.terms=TRUE)
#Partial effect plot for temperature indicates no effect on weight 

plot_predictions(mod1, condition=c("width", "year"))
plot_predictions(mod1, condition=list("width", "year"=2017:2018, "gear_temperature" = 1:1.99))

#Now make predictions while holding temperature at 1C
precollapse_dat %>%
  mutate(temp_bin = floor(gear_temperature)) %>%
  filter(temp_bin == 1) -> new.dat

p <- predict.gam(mod1, newdata=new.dat, type="response")
p <- as.data.frame(p)
pred <- cbind(new.dat, p)

#Not sure how this was done- there's no weight data at 1C in 2018 for predictions?

#what about adding julian day to attempt to correct for variation in molt timing?
mod1a <- gam(weight ~ s(width) + year + s(julian), data = precollapse_dat)
summary(mod1a) #negative coefficient for 2018, but fairly small effect
plot.gam(mod1a)
#plot julian day effect on response scale
plot_predictions(mod1a, condition = 'julian',
                 type = 'response', points = 0.5) +
  labs(y = "Expected response",
       title = "Average smooth effect of Julian day") +
  theme_bw()
#julian day effect can't be distinguished from zero


#now run gam for full dataset, including post collapse years 
mod2 <- gam(weight ~ s(width) + year + s(gear_temperature), data = full_dat)
summary(mod2) 
plot.gam(mod2, se=TRUE, shade=TRUE, all.terms=TRUE)

#and a model without year covariate, to fit a single regression across all years, 
  #and extract annual residuals 
mod3 <- gam(weight ~ s(width), data = full_dat)
summary(mod3) 

#extract annual L:W residuals from full model- I think this is the best way to
  #go about this?
full_dat$resid <- residuals(mod3)

full_dat %>% 
  group_by(year) %>%
  summarise(Avg_resid = mean(resid)) %>%
  ggplot(aes(year, Avg_resid)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(y = "Mean L:W residual") +
  ggtitle("Survey data: \nAnnual L:W residuals") 
#again, fairly small effect size in 2018

############################################################

#Now, let's calculate a relative condition factor to compare to fatty acid data for all
  #hepato sampled crab: Kn = observed weight/predicted weight (via regression fit from 
  #survey data, pooled across years) 

#Wrangle survey data - filtering out old shell and ovigerous females to match our 
  #collection criteria, and adding NBS data
ebs_haul %>%
  bind_rows(nbs_haul) %>% 
  rename_with(tolower) %>%
  mutate(year = as.numeric(str_extract(cruise, "\\d{4}"))) %>%
  filter(year %in% c(2011, 2015, 2017:2019, 2022, 2024),
         haul_type==3, 
         weight != "NA", 
         weight != 0,
         width != "NA",
         sex==1 & shell_condition<3 | sex==2 & shell_condition<3 & clutch_size==0) %>%
  mutate(log_weight = as.numeric(log(weight)),
         log_cw = as.numeric(log(width)), 
         year = as.factor(year)) -> dat2

#plot full dataset
dat2 %>%
  #filter(year %in% c(2022, 2024)) %>%
  ggplot(aes(width, weight,color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "gam", se = FALSE) 

#Predicted weight:fit L:W regression from 2011-2024, eliminating poorly sampled years
mod4 <- lm(log_weight ~ log_cw, data = dat2)
summary(mod4)
coef(mod4)
# log(W) = -8.213247  +  3.078925 * (L) on transformed scale    
# W = exp(-8.213247) * L^(3.078925)  on original scale 
#Note that we're overlooking a bias correction here should we need to back-transform!

#extract annual L:W residuals and plot
dat2$resid <- residuals(mod4)

dat2 %>% 
  group_by(year) %>%
  summarise(Avg_resid = mean(resid)) %>%
  ggplot(aes(year, Avg_resid)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(y = "Mean L:W residual") +
  ggtitle("Survey data: \nAnnual L:W residuals") 
#Very small residuals....

#now wrangle condition data for observed weights- no weights taken in 2024!
condition_master %>%
  filter(lme != "NA", #one crab collected outside the sampling design
         !vial_id %in% c("2019-65","2019-67","2019-68","2019-71","2019-66","2019-207", "2019-212"),#likely tanners
         maturity != 1,
         crab_wgt != "NA", 
        cw != "NA",) %>%
  mutate(log_weight = as.numeric(log(crab_wgt)),
         log_cw = as.numeric(log(cw)), 
         year = as.factor(year))-> cond_dat

#predict weight of condition crab using regression fit with survey weight at size data
pred_dat <-  cond_dat %>%
  dplyr::select(log_cw, log_weight, year)

predictions <- predict(mod4, pred_dat, interval = 'confidence') %>% 
  as.data.frame()

bind_cols(cond_dat, predictions) %>%
  mutate(k_n = log_weight/fit,
         resid = log_weight - fit) -> dat_check

dat_check %>% 
  group_by(year) %>%
  summarise(Avg_resid = mean(resid)) %>%
  ggplot(aes(year, Avg_resid)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(y = "Mean L:W residual") +
  ggtitle("Condition crab: \nAnnual L:W residuals") 


#plot relative condition, k_n, vrs total FA
dat_check %>%
ggplot(aes(k_n, Total_FA_Conc_WWT, color=factor(year))) +
  geom_point() +
  theme_bw()  +
  #geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Relative Condition (observed/predicted weight at size)", y = "Total FA per WWT (mg FA/g WWT)") +
  ggtitle("Study samples: \nFA vrs Morphometric Condition") 

#model relationship between k_n and total FA
mod5 <- lm(k_n ~ Total_FA_Conc_WWT + year, data = dat_check)
summary(mod5)
#pretty darn lousy....

#now plot L:W residuals vrs total FA
dat_check %>%
  ggplot(aes(resid, Total_FA_Conc_WWT, color=factor(year))) +
  geom_point() +
  theme_bw()  +
  #geom_smooth(method = "lm", se = FALSE) +
  labs(x= "L:W Residual", y = "Total FA per WWT (mg FA/g WWT)") +
  ggtitle("Study samples: \nA Condition vrs L:W residuals") 

#Fultons K Condition factor vrs % WWT
cond_dat %>%
  filter(!vial_id %in% c("2023-147", "2022-AKK-175")) %>% #outliers based on wgt- likely back deck errors
  mutate(K=crab_wgt/(cw^3)) %>%
  ggplot(aes(K, Total_FA_Conc_WWT, color=factor(year))) +
  geom_point() +
  theme_bw()  +
  #geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Fultons K Condition Factor", y = "Total FA per DWT (mg FA/g WWT)") 
#very little variation in weight at size as compared to total FA, but Fultons K
#doesn't account for variation in crab size well

#Follow up: males vrs females modeled separately (2018 sample size discrepancy
  #due to filtering out females maybe?)








  
  
  
  
  














