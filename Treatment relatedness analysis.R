# This code was used to investigate the effect of treatment on focal species reproductive output (Q1)

# load relevant libraries
library(ggpubr)
library(tidyr)
library(dplyr)
library(readr)
library(gridExtra)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(ggstance)
library(MuMIn)
library(visreg)
library(MASS)
library(jtools)
library(interactions)
library(broom.mixed)
library(devtools)
library(lme4)
library(nlme)
library(sm)
library(performance)
library(dfoptim)
library(optimx)
library(multcomp)
library(AICcmodavg)
library(ggExtra)
library(sjPlot)

# import data
experimental_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/experimental_data_plot_summary.csv")

## remove unwanted values from the dataframe
# factorise treatment variables for experimental_dat
experimental_dat$treatment<-as.factor(experimental_dat$treatment)
# keep experimental_dat only with correct treatments
experimental_dat<-experimental_dat[which(experimental_dat$correct_treatment == "Y"), ]
# remove problematic experimental_dat points (plots) as verified with testing
experimental_dat<-experimental_dat[which(experimental_dat$problematic == "N"), ]
# create df for plants dead at final seed collection/count
experimental_dat<-experimental_dat[which(experimental_dat$alive_at_collection == "N"), ]

# create seperate df for each species (all plots)
arca_exp_dat<-experimental_dat[which(experimental_dat$species == "arca"), ]
momo_exp_dat<-experimental_dat[which(experimental_dat$species == "momo"), ]
peai_exp_dat<-experimental_dat[which(experimental_dat$species == "peai"), ]
avba_exp_dat<-experimental_dat[which(experimental_dat$species == "avba"), ]

#### Note here that previous investigation found negative binomial to be a better fit than poisson (in which dispersion was detected), this was the case for all species (see R code for "GLM distribution testing")

## ARCA

# set the reference treatment level
arca_exp_dat$treatment <- relevel(arca_exp_dat$treatment,"control")
arca_exp_dat$treatment <- relevel(arca_exp_dat$treatment,"thinned_control")

# plot the treatment data as a response to number of inflorescences
ggplot()+
  geom_boxplot(data=arca_exp_dat, aes(x=treatment, y=no_inflor_with_seeds_total))

## test only fixed effects
arca_fixed<-glm.nb(no_inflor_with_seeds_total~treatment,
                   data=arca_exp_dat)
summary(arca_fixed)
AICc(arca_fixed) #357.8405

## test random slopes
# number of competitors within plot
arca_comp<-glmer.nb(no_inflor_with_seeds_total~treatment+(no_comps_total|plot),
                      data=arca_exp_dat)
summary(arca_comp)
AICc(arca_comp) #361.1171

# richness within plot
arca_richness<-glmer.nb(no_inflor_with_seeds_total~treatment+(richness|plot),
                      data=arca_exp_dat)
summary(arca_richness)
AICc(arca_richness) #357.6867

# plot within block
arca_plotwblock<-glmer.nb(no_inflor_with_seeds_total~treatment+(plot|block),
                      data=arca_exp_dat)
summary(arca_plotwblock)
AICc(arca_plotwblock) #362.6783

# presence of perrenials within plot
arca_perrenials<-glmer.nb(no_inflor_with_seeds_total~treatment+(correct_ID|plot),
                 data=arca_exp_dat)
summary(arca_perrenials)
AICc(arca_perrenials) #361.0439

## test random intercepts
# block
arca_block<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|block),
               data=arca_exp_dat)
summary(arca_block)
AICc(arca_block) #357.263

# plot
arca_plot<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|plot),
                     data=arca_exp_dat)
summary(arca_plot)
AICc(arca_plot) #355.9376

# subplot
arca_subplot<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|subplot),
                    data=arca_exp_dat)
summary(arca_subplot)
AICc(arca_subplot) #360.6262

# plot within block + block
arca_block_plot<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|block/plot),
                       data=arca_exp_dat)
summary(arca_block_plot)
AICc(arca_block_plot) #358.5431

# subplot within plot + plot
arca_plot_subplot<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|plot/subplot),
                          data=arca_exp_dat)
summary(arca_plot_subplot)
AICc(arca_plot_subplot) #358.8265

# final arca model:
arca_final<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|plot),
                                data=arca_exp_dat)
summary(arca_final)
AICc(arca_final) #355.9376

# use visreg function to visualise the effect of treatment on performance based on model output 
visreg(arca_final, "treatment", scale='response', rug=FALSE)

## MOMO

# set the reference treatment level
momo_exp_dat$treatment <- relevel(momo_exp_dat$treatment,"control")
momo_exp_dat$treatment <- relevel(momo_exp_dat$treatment,"thinned_control")

# plot the treatment data as a response to number of seeds
ggplot()+
  geom_boxplot(data=momo_exp_dat, aes(x=treatment, y=developed_seeds_counted))

## test only fixed effects
momo_fixed<-glm.nb(developed_seeds_counted~treatment,
                   data=momo_exp_dat)
summary(momo_fixed)
AICc(momo_fixed) #535.3429

## test random slopes
# number of competitors within plot
momo_comp<-glmer.nb(developed_seeds_counted~treatment+(no_comps_total|plot),
                    data=momo_exp_dat)
summary(momo_comp)
AICc(momo_comp) #533.5212

# richness within plot
momo_richness<-glmer.nb(developed_seeds_counted~treatment+(richness|plot),
                        data=momo_exp_dat)
summary(momo_richness)
AICc(momo_richness) #532.5725

# plot within block
momo_plotwblock<-glmer.nb(developed_seeds_counted~treatment+(plot|block),
                          data=momo_exp_dat)
summary(momo_plotwblock)
AICc(momo_plotwblock) #536.8556

# presence of perrenials within plot
momo_perrenials<-glmer.nb(developed_seeds_counted~treatment+(correct_ID|plot),
                          data=momo_exp_dat)
summary(momo_perrenials)
AICc(momo_perrenials) #533.8536

## test random intercepts
# block
momo_block<-glmer.nb(developed_seeds_counted~treatment+(1|block),
                     data=momo_exp_dat)
summary(momo_block)
AICc(momo_block) #536.0341

# plot
momo_plot<-glmer.nb(developed_seeds_counted~treatment+(1|plot),
                    data=momo_exp_dat)
summary(momo_plot)
AICc(momo_plot) #528.3982

# subplot
momo_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|subplot),
                       data=momo_exp_dat)
summary(momo_subplot)
AICc(momo_subplot) #538.0852

# plot within block + block
momo_block_plot<-glmer.nb(developed_seeds_counted~treatment+(1|block/plot),
                          data=momo_exp_dat)
summary(momo_block_plot)
AICc(momo_block_plot) #531.2367

# subplot within plot + plot
momo_plot_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|plot/subplot),
                            data=momo_exp_dat)
summary(momo_plot_subplot)
AICc(momo_plot_subplot) #531.2367

# final momo model:
momo_final<-glmer.nb(developed_seeds_counted~treatment+(1|plot),
                     data=momo_exp_dat)
summary(momo_final)
AICc(momo_final) #528.3982

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(momo_final, "treatment", scale='response', rug=FALSE)

## PEAI

# set the reference treatment level
peai_exp_dat$treatment <- relevel(peai_exp_dat$treatment,"control")
peai_exp_dat$treatment <- relevel(peai_exp_dat$treatment,"thinned_control")

# plot the treatment data as a response to number of seeds
ggplot()+
  geom_boxplot(data=peai_exp_dat, aes(x=treatment, y=developed_seeds_counted))

## test only fixed effects
peai_fixed<-glm.nb(developed_seeds_counted~treatment,
                   data=peai_exp_dat)
summary(peai_fixed)
AICc(peai_fixed) #831.1842

## test random slopes
# number of competitors within plot
peai_comp<-glmer.nb(developed_seeds_counted~treatment+(no_comps_total|plot),
                    data=peai_exp_dat)
summary(peai_comp)
AICc(peai_comp) #839.9845

# richness within plot
peai_richness<-glmer.nb(developed_seeds_counted~treatment+(richness|plot),
                        data=peai_exp_dat)
summary(peai_richness)
AICc(peai_richness) #839.9845

# plot within block
peai_plotwblock<-glmer.nb(developed_seeds_counted~treatment+(plot|block),
                          data=peai_exp_dat)
summary(peai_plotwblock)
AICc(peai_plotwblock) #839.9845

# presence of perrenials within plot
peai_perrenials<-glmer.nb(developed_seeds_counted~treatment+(correct_ID|plot),
                          data=peai_exp_dat)
summary(peai_perrenials)
AICc(peai_perrenials) # did not converge

## test random intercepts
# block
peai_block<-glmer.nb(developed_seeds_counted~treatment+(1|block),
                     data=peai_exp_dat)
summary(peai_block)
AICc(peai_block) #833.9954

# plot
peai_plot<-glmer.nb(developed_seeds_counted~treatment+(1|plot),
                    data=peai_exp_dat)
summary(peai_plot)
AICc(peai_plot) #833.9954

# subplot
peai_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|subplot),
                       data=peai_exp_dat)
summary(peai_subplot)
AICc(peai_subplot) #833.9954

# plot within block + block
peai_block_plot<-glmer.nb(developed_seeds_counted~treatment+(1|block/plot),
                          data=peai_exp_dat)
summary(peai_block_plot)
AICc(peai_block_plot) #836.9262

# subplot within plot + plot
peai_plot_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|plot/subplot),
                            data=peai_exp_dat)
summary(peai_plot_subplot)
AICc(peai_plot_subplot) #836.9262

# final peai model:
peai_final<-glm.nb(developed_seeds_counted~treatment,
                   data=peai_exp_dat)
summary(peai_final)
AICc(peai_final) #831.1842

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(peai_final, "treatment", scale='response', rug=FALSE)

## AVBA

# set the reference treatment level
avba_exp_dat$treatment <- relevel(avba_exp_dat$treatment,"control")
avba_exp_dat$treatment <- relevel(avba_exp_dat$treatment,"thinned_control")

# plot the treatment data as a response to number of seeds
ggplot()+
  geom_boxplot(data=avba_exp_dat, aes(x=treatment, y=developed_seeds_counted))

## test only fixed effects
avba_fixed<-glm.nb(developed_seeds_counted~treatment,
                   data=avba_exp_dat)
summary(avba_fixed)
AICc(avba_fixed) #513.8866

## test random slopes
# number of competitors within plot
avba_comp<-glmer.nb(developed_seeds_counted~treatment+(no_comps_total|plot),
                    data=avba_exp_dat)
summary(avba_comp)
AICc(avba_comp) #517.9078

# richness within plot
avba_richness<-glmer.nb(developed_seeds_counted~treatment+(richness|plot),
                        data=avba_exp_dat)
summary(avba_richness)
AICc(avba_richness) #522.2536

# plot within block
avba_plotwblock<-glmer.nb(developed_seeds_counted~treatment+(plot|block),
                          data=avba_exp_dat)
summary(avba_plotwblock)
AICc(avba_plotwblock) #522.4496

# presence of perrenials within plot
avba_perrenials<-glmer.nb(developed_seeds_counted~treatment+(correct_ID|plot),
                          data=avba_exp_dat)
summary(avba_perrenials)
AICc(avba_perrenials) #521.7081

## test random intercepts
# block
avba_block<-glmer.nb(developed_seeds_counted~treatment+(1|block),
                     data=avba_exp_dat)
summary(avba_block)
AICc(avba_block) #516.6311

# plot
avba_plot<-glmer.nb(developed_seeds_counted~treatment+(1|plot),
                    data=avba_exp_dat)
summary(avba_plot)
AICc(avba_plot) #516.6311

# subplot
avba_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|subplot),
                       data=avba_exp_dat)
summary(avba_subplot)
AICc(avba_subplot) #516.5862

# plot within block + block
avba_block_plot<-glmer.nb(developed_seeds_counted~treatment+(1|block/plot),
                          data=avba_exp_dat)
summary(avba_block_plot)
AICc(avba_block_plot) #519.4833

# subplot within plot + plot
avba_plot_subplot<-glmer.nb(developed_seeds_counted~treatment+(1|plot/subplot),
                            data=avba_exp_dat)
summary(avba_plot_subplot)
AICc(avba_plot_subplot) #519.4833

# final avba model:
avba_final<-glm.nb(developed_seeds_counted~treatment,
                   data=avba_exp_dat)
summary(avba_final)
AICc(avba_final) #513.8866

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(avba_final, "treatment", scale='response', rug=FALSE)

################### FINAL MODELS

## ARCA
# relvel treatments
arca_exp_dat$treatment <- relevel(arca_exp_dat$treatment,"control")
arca_exp_dat$treatment <- relevel(arca_exp_dat$treatment,"thinned_control")
# final model
arca_final<-glmer.nb(no_inflor_with_seeds_total~treatment+(1|plot),
                     data=arca_exp_dat)
summary(arca_final)

## MOMO
# relvel treatments
momo_exp_dat$treatment <- relevel(momo_exp_dat$treatment,"control")
momo_exp_dat$treatment <- relevel(momo_exp_dat$treatment,"thinned_control")
# final model
momo_final<-glmer.nb(developed_seeds_counted~treatment+(1|plot),
                     data=momo_exp_dat)
summary(momo_final)

## PEAI
# relvel treatments
peai_exp_dat$treatment <- relevel(peai_exp_dat$treatment,"control")
peai_exp_dat$treatment <- relevel(peai_exp_dat$treatment,"thinned_control")
# final model
peai_final<-glm.nb(developed_seeds_counted~treatment,
                   data=peai_exp_dat)
summary(peai_final)

## AVBA
# relvel treatments
avba_exp_dat$treatment <- relevel(avba_exp_dat$treatment,"control")
avba_exp_dat$treatment <- relevel(avba_exp_dat$treatment,"thinned_control")
# final model
avba_final<-glm.nb(developed_seeds_counted~treatment,
                   data=avba_exp_dat)
summary(avba_final)
