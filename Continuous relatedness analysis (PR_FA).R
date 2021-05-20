# This code was used to investigate the effect of ANPD (calculated including abundance) on focal species reproductive output (Q1)

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
# create df for plants dead at final seed collection/count
experimental_dat<-experimental_dat[which(experimental_dat$alive_at_collection == "N"), ]
#create df excluding monoculture and no competitor treatment
experimental_dat<-experimental_dat[-which(experimental_dat$treatment == "no_competitors" | experimental_dat$treatment == "monoculture"), ]

# create seperate df for each species (all plots)
arca_exp_dat<-experimental_dat[which(experimental_dat$species == "arca"), ]
momo_exp_dat<-experimental_dat[which(experimental_dat$species == "momo"), ]
peai_exp_dat<-experimental_dat[which(experimental_dat$species == "peai"), ]
avba_exp_dat<-experimental_dat[which(experimental_dat$species == "avba"), ]

#### Note here that previous investigation found negative binomial to be a better fit than poisson (in which dispersion was detected), this was the case for all species (see R code for "GLM distribution testing")

## ARCA

## test only fixed effects
arca_fixed<-glm.nb(no_inflor_with_seeds_total~PR_FA,
                   data=arca_exp_dat)
summary(arca_fixed)
AICc(arca_fixed) #312.4002

## test random slopes
# number of competitors within plot
arca_comp<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(no_comps_total|plot),
                    data=arca_exp_dat)
summary(arca_comp)
AICc(arca_comp) #315.6693

# richness within plot
arca_richness<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(richness|plot),
                        data=arca_exp_dat)
summary(arca_richness)
AICc(arca_richness) #313.3471

# plot within block
arca_plotwblock<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(plot|block),
                          data=arca_exp_dat)
summary(arca_plotwblock)
AICc(arca_plotwblock) #317.4286

# presence of perrenials within plot
arca_perrenials<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(correct_ID|plot),
                          data=arca_exp_dat)
summary(arca_perrenials)
AICc(arca_perrenials) #315.1856

## test random intercepts
# block
arca_block<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|block),
                     data=arca_exp_dat)
summary(arca_block)
AICc(arca_block) #312.7103

# plot
arca_plot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|plot),
                    data=arca_exp_dat)
summary(arca_plot)
AICc(arca_plot) #311.0124

# subplot
arca_subplot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|subplot),
                       data=arca_exp_dat)
summary(arca_subplot)
AICc(arca_subplot) #314.7166

# plot within block + block
arca_block_plot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|block/plot),
                          data=arca_exp_dat)
summary(arca_block_plot)
AICc(arca_block_plot) #313.4281

# subplot within plot + plot
arca_plot_subplot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|plot/subplot),
                            data=arca_exp_dat)
summary(arca_plot_subplot)
AICc(arca_plot_subplot) #313.4197

# final arca model:
arca_final<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=arca_exp_dat)
summary(arca_final)
AICc(arca_final) #311.0124

# use visreg function to visualise the effect of PR_FA on performance based on model output 
visreg(arca_final, "PR_FA", scale='response', rug=FALSE)

## MOMO

## test only fixed effects
momo_fixed<-glm.nb(developed_seeds_counted~PR_FA,
                   data=momo_exp_dat)
summary(momo_fixed)
AICc(momo_fixed) #458.5519

## test random slopes
# number of competitors within plot
momo_comp<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(no_comps_total|plot),
                    data=momo_exp_dat)
summary(momo_comp)
AICc(momo_comp) #457.6076

# richness within plot
momo_richness<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(richness|plot),
                        control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                        data=momo_exp_dat)
summary(momo_richness)
AICc(momo_richness) #456.7671

# plot within block
momo_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(plot|block),
                          data=momo_exp_dat)
summary(momo_plotwblock)
AICc(momo_plotwblock) #462.9597

# presence of perrenials within plot
momo_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(correct_ID|plot),
                          data=momo_exp_dat)
summary(momo_perrenials)
AICc(momo_perrenials) #457.8295

## test random intercepts
# block
momo_block<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block),
                     data=momo_exp_dat)
summary(momo_block)
AICc(momo_block) #460.047

# plot
momo_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot),
                    control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                    data=momo_exp_dat)
summary(momo_plot)
AICc(momo_plot) #453.1377

# subplot
momo_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|subplot),
                       data=momo_exp_dat)
summary(momo_subplot)
AICc(momo_subplot) #460.8563

# plot within block + block
momo_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block/plot),
                          data=momo_exp_dat)
summary(momo_block_plot)
AICc(momo_block_plot) #455.528

# subplot within plot + plot
momo_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot/subplot),
                            control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                            data=momo_exp_dat)
summary(momo_plot_subplot)
AICc(momo_plot_subplot) #454.988

# final momo model:
momo_final<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=momo_exp_dat)
summary(momo_final)
AICc(momo_final) #453.1377

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(momo_final, "PR_FA", scale='response', rug=FALSE)

## PEAI

## test only fixed effects
peai_fixed<-glm.nb(developed_seeds_counted~PR_FA,
                   data=peai_exp_dat)
summary(peai_fixed)
AICc(peai_fixed) #637.5767

## test random slopes
# number of competitors within plot
peai_comp<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(no_comps_total|plot),
                    data=peai_exp_dat)
summary(peai_comp)
AICc(peai_comp) #645.2019

# richness within plot
peai_richness<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(richness|plot),
                        data=peai_exp_dat)
summary(peai_richness)
AICc(peai_richness) #645.2019

# plot within block
peai_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(plot|block),
                          data=peai_exp_dat)
summary(peai_plotwblock)
AICc(peai_plotwblock) #645.2019

# presence of perrenials within plot
peai_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(correct_ID|plot),
                          data=peai_exp_dat)
summary(peai_perrenials)
AICc(peai_perrenials) # did not converge

## test random intercepts
# block
peai_block<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block),
                     data=peai_exp_dat)
summary(peai_block)
AICc(peai_block) #639.9914

# plot
peai_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot),
                    data=peai_exp_dat)
summary(peai_plot)
AICc(peai_plot) #639.9914

# subplot
peai_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|subplot),
                       data=peai_exp_dat)
summary(peai_subplot)
AICc(peai_subplot) #639.9914

# plot within block + block
peai_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block/plot),
                          data=peai_exp_dat)
summary(peai_block_plot)
AICc(peai_block_plot) #642.5298

# subplot within plot + plot
peai_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot/subplot),
                            data=peai_exp_dat)
summary(peai_plot_subplot)
AICc(peai_plot_subplot) #642.5298

# final peai model:
peai_final<-glm.nb(developed_seeds_counted~PR_FA,
                   data=peai_exp_dat)
summary(peai_final)
AICc(peai_final) #831.1842

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(peai_final, "PR_FA", scale='response', rug=FALSE)

## AVBA

## test only fixed effects
avba_fixed<-glm.nb(developed_seeds_counted~scale(PR_FA),
                   data=avba_exp_dat)
summary(avba_fixed)
AICc(avba_fixed) #388.9743

## test random slopes
# number of competitors within plot
avba_comp<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(no_comps_total|plot),
                    data=avba_exp_dat)
summary(avba_comp)
AICc(avba_comp) #382.6051

# richness within plot
avba_richness<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(richness|plot),
                        data=avba_exp_dat)
summary(avba_richness)
AICc(avba_richness) #393.5177

# plot within block
avba_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(plot|block),
                          data=avba_exp_dat)
summary(avba_plotwblock)
AICc(avba_plotwblock) #395.8264

# presence of perrenials within plot
avba_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(correct_ID|plot),
                          data=avba_exp_dat)
summary(avba_perrenials)
AICc(avba_perrenials) #396.0805 (failed to converge)

## test random intercepts
# block
avba_block<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block),
                     data=avba_exp_dat)
summary(avba_block)
AICc(avba_block) #391.3686

# plot
avba_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot),
                    data=avba_exp_dat)
summary(avba_plot)
AICc(avba_plot) #391.3249

# subplot
avba_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|subplot),
                       data=avba_exp_dat)
summary(avba_subplot)
AICc(avba_subplot) #391.3686

# plot within block + block
avba_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|block/plot),
                          data=avba_exp_dat)
summary(avba_block_plot)
AICc(avba_block_plot) #393.836

# subplot within plot + plot
avba_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot/subplot),
                            data=avba_exp_dat)
summary(avba_plot_subplot)
AICc(avba_plot_subplot) #393.836

# final avba model:
avba_final<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(no_comps_total|plot),
                     data=avba_exp_dat)
summary(avba_final)
AICc(avba_final) #374.7471

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(avba_final, "PR_FA", scale='response', rug=FALSE)

################### FINAL MODELS

# ARCA
arca_final<-glmer.nb(no_inflor_with_seeds_total~scale(PR_FA)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=arca_exp_dat)
summary(arca_final)

# MOMO
momo_final<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=momo_exp_dat)
summary(momo_final)

# PEAI
peai_final<-glm.nb(developed_seeds_counted~scale(PR_FA),
                   data=peai_exp_dat)
summary(peai_final)

# AVBA
avba_final<-glmer.nb(developed_seeds_counted~scale(PR_FA)+(no_comps_total|plot),
                   data=avba_exp_dat)
summary(avba_final)

