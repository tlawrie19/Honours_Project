# This code was used to investigate which distributional family (for glm) was appropriate for each species (either poisson and glm)

# load relevant libraries
library(devtools)
library(lme4)
library(nlme)
library(ggplot2)  
library(sm)
library(performance)
library(MASS)
library(dfoptim)
library(optimx)
library(multcomp)
library(visreg)
library(AICcmodavg)

# import data
experimental_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/experimental_data_plot_summary.csv")

## remove unwanted values from the dataframe
# factorise treatment variables for experimental_dat
experimental_dat$treatment<-as.factor(experimental_dat$treatment)
# keep experimental_dat only with correct treatments
experimental_dat<-experimental_dat[which(experimental_dat$correct_treatment == "Y"), ]
# create df for plants dead at final seed collection/count
experimental_dat<-experimental_dat[which(experimental_dat$alive_at_collection == "N"), ]
#create df excluding monoculture and no competitor treatment
experimental_dat<-experimental_dat[-which(experimental_dat$treatment == "no_competitors" | experimental_dat$treatment == "monoculture"), ]

# create seperate df for each species (all plots)
arca_exp_dat<-experimental_dat[which(experimental_dat$species == "arca"), ]
momo_exp_dat<-experimental_dat[which(experimental_dat$species == "momo"), ]
peai_exp_dat<-experimental_dat[which(experimental_dat$species == "peai"), ]
avba_exp_dat<-experimental_dat[which(experimental_dat$species == "avba"), ]

################################### ARCA

# factorise treatment variables
arca_exp_dat$treatment<-as.factor(arca_exp_dat$treatment)

## assesing distributional family requirements (with glm)
# poisson
arca_poisson_glm <- glm(no_inflor_with_seeds_total~treatment, data=arca_exp_dat,family=poisson)
print(summary(arca_poisson_glm)) ### AIC = 371.84
AICc(arca_poisson_glm) # 374.4597

# negative binomial
arca_nb_glm <- glm.nb(no_inflor_with_seeds_total~treatment, data=arca_exp_dat)
print(summary(arca_nb_glm)) ### AIC = 338.19
AICc(arca_nb_glm) # 341.42

# check overdispersion of poisson
check_overdispersion(arca_poisson_glm) #overdispersion is detected

# plot residuals
plot(arca_poisson_glm) 
plot(arca_nb_glm) #### nb better AIC and residual variance 

################################### MOMO

# factorise treatment variables
momo_exp_dat$treatment<-as.factor(momo_exp_dat$treatment)

## assesing distributional family requirements (with glm)
# poisson
momo_poisson_glm <- glm(developed_seeds_counted~treatment, data=momo_exp_dat,family=poisson)
print(summary(momo_poisson_glm)) ### AIC = 1544
AICc(momo_poisson_glm)

# negative binomial
momo_nb_glm <- glm.nb(developed_seeds_counted~treatment, data=momo_exp_dat)
print(summary(momo_nb_glm)) ### AIC = 519.68
AICc(momo_nb_glm)

# check overdispersion of poisson
check_overdispersion(momo_poisson_glm) #overdispersion is detected

# plot residuals
plot(momo_poisson_glm) 
plot(momo_nb_glm) #### nb better AIC and residual variance 

################################### PEAI

# factorise treatment variables
peai_exp_dat$treatment<-as.factor(peai_exp_dat$treatment)

## assesing distributional family requirements (with glm)
# poisson
peai_poisson_glm <- glm(developed_seeds_counted~treatment, data=peai_exp_dat,family=poisson)
print(summary(peai_poisson_glm)) ### AIC = 1544
AICc(peai_poisson_glm)

# negative binomial
peai_nb_glm <- glm.nb(developed_seeds_counted~treatment, data=peai_exp_dat)
print(summary(peai_nb_glm)) ### AIC = 519.68
AICc(peai_nb_glm)

# check overdispersion of poisson
check_overdispersion(peai_poisson_glm) #overdispersion is detected

# plot residuals
plot(peai_poisson_glm) 
plot(peai_nb_glm) #### nb better AIC and residual variance 

################################### AVBA

# factorise treatment variables
avba_exp_dat$treatment<-as.factor(avba_exp_dat$treatment)

## assesing distributional family requirements (with glm)
# poisson
avba_poisson_glm <- glm(developed_seeds_counted~treatment, data=avba_exp_dat,family=poisson)
print(summary(avba_poisson_glm)) ### AIC = 1544
AICc(avba_poisson_glm)

# negative binomial
avba_nb_glm <- glm.nb(developed_seeds_counted~treatment, data=avba_exp_dat)
print(summary(avba_nb_glm)) ### AIC = 519.68
AICc(avba_nb_glm)

# check overdispersion of poisson
check_overdispersion(avba_poisson_glm) #overdispersion is detected

# plot residuals
plot(avba_poisson_glm) 
plot(avba_nb_glm) #### nb better AIC and residual variance 
