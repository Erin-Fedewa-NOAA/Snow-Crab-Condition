#Goal: Assess the relationship b/w total fatty acid concentration and indirect proxies for condition

#NOTE: We're removing 2019 from validation models b/c methods differed (i.e. hepatos were dissected
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
source("./script/stan_utils.R")

#load data 
condition_master <- read.csv("./data/total_FA_master.csv")

#plotting
my_colors <- RColorBrewer::brewer.pal(7, "GnBu")[c(3,5,7)]
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

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
            family = student, #Student's t distribution- more robust to outliers
            Total_FA_Conc_DWT ~ 1 + Perc_DWT,
            seed=1,
            save_pars = save_pars(all = TRUE),
            control = list(adapt_delta = 0.999, max_treedepth = 14))

#Save output
saveRDS(condition_1, file = "./output/condition_1.rds")
condition_1 <- readRDS("./output/condition_1.rds")

summary(condition_1)
bayes_R2(condition_1) #.67
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

#Plot model fit - Fig 2 for Manuscript 
new.dat %>%
  ggplot(aes(x = Perc_DWT, y = Total_FA_Conc_DWT)) +
  #geom_abline(intercept = fixef(condition_1)[1], 
             # slope     = fixef(condition_1)[2],
             # size = .8, color = "black") +
  geom_smooth(data = mu, aes(y=Estimate, ymin= Q2.5, ymax= Q97.5), 
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
ggsave("./figures/Fig2.png", height = 4, width = 5, units = "in", dpi = 300)


############################
#Morphometric condition metrics

#Crab weight at size vrs % DWT by sex
new.dat %>%
  filter(!vial_id %in% c("2023-147", "2022-AKK-175")) %>% #outliers based on wgt- likely back deck errors
  mutate(lw = crab_wgt/cw) %>%
  ggplot(aes(lw, Total_FA_Conc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Crab weight/size ratio", y = "Total FA per DWT (mg FA/g DWT)")

#% LM: weight/size ratio vrs total FA concentration - no 2019
new.dat %>%
  filter(year > 2019,
         !vial_id %in% c("2023-147", "2022-AKK-175")) %>% #outliers based on wgt- likely back deck errors-> new.dat
  mutate(lw = crab_wgt/cw) -> wgt_dat
mod2 <- lm(Total_FA_Conc_DWT~lw, data=wgt_dat) 
summary(mod2)

# Condition factor K vrs % DWT
new.dat %>%
  filter(!vial_id %in% c("2023-147", "2022-AKK-175")) %>% #outliers based on wgt- likely back deck errors
  mutate(K=crab_wgt/(cw^3)) %>%
  ggplot(aes(K, Total_FA_Conc_DWT, color=factor(year))) +
  geom_point() +
  theme_bw()  +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x= "Fultons K Condition Factor", y = "Total FA per DWT (mg FA/g WWT)") 