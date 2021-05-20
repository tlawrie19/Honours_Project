# This code was used to investigate the effect of ANPD (calculated excluded abundance) on focal species reproductive output (Q1)
# Further, this code was used to test for an association between native range ANPD values and fitness outcomes from the invaded range (Q3)

# load libraries
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
library(emmeans)

# import data
native_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Native Range Data/native_range_plot_summary.csv")
experimental_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/experimental_data_plot_summary.csv")
pj_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Native Range Data/ANPD_native_v_invaded_summary.csv")

## create seperate df for each species native range data (with condition + plots only <100)
arca_native_dat<-native_dat[which(native_dat$focal_species=="arca"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]
momo_native_dat<-native_dat[which(native_dat$focal_species=="momo"&native_dat$condition=="H"&native_dat$plot_size_m2<=100&native_dat$plot_size_m2>0),]
peai_native_dat<-native_dat[which(native_dat$focal_species=="peai"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]
avba_native_dat<-native_dat[which(native_dat$focal_species=="avba"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]

## create seperate df for each species PJ data
arca_pj_dat<-pj_dat[which(pj_dat$focal_spec=="arca"),]
momo_pj_dat<-pj_dat[which(pj_dat$focal_spec=="momo"),]
peai_pj_dat<-pj_dat[which(pj_dat$focal_spec=="peai"),]
avba_pj_dat<-pj_dat[which(pj_dat$focal_spec=="avba"),]


### Experimental data continuous phylogenetic values

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
arca_fixed<-glm.nb(no_inflor_with_seeds_total~PR_F,
                   data=arca_exp_dat)
summary(arca_fixed)
AICc(arca_fixed) #311.9634

## test random slopes
# number of competitors within plot
arca_comp<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(no_comps_total|plot),
                    data=arca_exp_dat)
summary(arca_comp)
AICc(arca_comp) #315.4746

# richness within plot
arca_richness<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(richness|plot),
                        data=arca_exp_dat)
summary(arca_richness)
AICc(arca_richness) #313.368

# plot within block
arca_plotwblock<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(plot|block),
                          data=arca_exp_dat)
summary(arca_plotwblock)
AICc(arca_plotwblock) #317.0456

# presence of perrenials within plot
arca_perrenials<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(correct_ID|plot),
                          data=arca_exp_dat)
summary(arca_perrenials)
AICc(arca_perrenials) #314.9411

## test random intercepts
# block
arca_block<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|block),
                     data=arca_exp_dat)
summary(arca_block)
AICc(arca_block) #312.3714

# plot
arca_plot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|plot),
                    data=arca_exp_dat)
summary(arca_plot)
AICc(arca_plot) #310.7337

# subplot
arca_subplot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|subplot),
                       data=arca_exp_dat)
summary(arca_subplot)
AICc(arca_subplot) #314.2798

# plot within block + block
arca_block_plot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|block/plot),
                          data=arca_exp_dat)
summary(arca_block_plot)
AICc(arca_block_plot) #313.1409

# subplot within plot + plot
arca_plot_subplot<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|plot/subplot),
                            data=arca_exp_dat)
summary(arca_plot_subplot)
AICc(arca_plot_subplot) #313.1409

# final arca model:
arca_final<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|plot),
                     data=arca_exp_dat)
summary(arca_final)
AICc(arca_final) #310.7337

# use visreg function to visualise the effect of PR_F on performance based on model output 
visreg(arca_final, "PR_F", scale='response', rug=FALSE)

## MOMO

## test only fixed effects
momo_fixed<-glm.nb(developed_seeds_counted~PR_F,
                   data=momo_exp_dat)
summary(momo_fixed)
AICc(momo_fixed) #458.6371

## test random slopes
# number of competitors within plot
momo_comp<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                    data=momo_exp_dat)
summary(momo_comp)
AICc(momo_comp) #457.5703

# richness within plot
momo_richness<-glmer.nb(developed_seeds_counted~scale(PR_F)+(richness|plot),
                        control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                        data=momo_exp_dat)
summary(momo_richness)
AICc(momo_richness) #456.7774

# plot within block
momo_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_F)+(plot|block),
                          data=momo_exp_dat)
summary(momo_plotwblock)
AICc(momo_plotwblock) #463.0471

# presence of perrenials within plot
momo_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_F)+(correct_ID|plot),
                          data=momo_exp_dat)
summary(momo_perrenials)
AICc(momo_perrenials) #457.7744

## test random intercepts
# block
momo_block<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block),
                     data=momo_exp_dat)
summary(momo_block)
AICc(momo_block) #460.1059

# plot
momo_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                    control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                    data=momo_exp_dat)
summary(momo_plot)
AICc(momo_plot) #453.0925

# subplot
momo_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|subplot),
                       data=momo_exp_dat)
summary(momo_subplot)
AICc(momo_subplot) #460.9414

# plot within block + block
momo_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block/plot),
                          data=momo_exp_dat)
summary(momo_block_plot)
AICc(momo_block_plot) #455.4832

# subplot within plot + plot
momo_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot/subplot),
                            control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                            data=momo_exp_dat)
summary(momo_plot_subplot)
AICc(momo_plot_subplot) #454.9842

# final momo model:
momo_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=momo_exp_dat)
summary(momo_final)
AICc(momo_final) #453.0925

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(momo_final, "PR_F", scale='response', rug=FALSE)

## PEAI

## test only fixed effects
peai_fixed<-glm.nb(developed_seeds_counted~PR_F,
                   data=peai_exp_dat)
summary(peai_fixed)
AICc(peai_fixed) #637.2174

## test random slopes
# number of competitors within plot
peai_comp<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                    data=peai_exp_dat)
summary(peai_comp)
AICc(peai_comp) #644.8425

# richness within plot
peai_richness<-glmer.nb(developed_seeds_counted~scale(PR_F)+(richness|plot),
                        data=peai_exp_dat)
summary(peai_richness)
AICc(peai_richness) #644.8425

# plot within block
peai_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_F)+(plot|block),
                          data=peai_exp_dat)
summary(peai_plotwblock)
AICc(peai_plotwblock) #644.8425

# presence of perrenials within plot
peai_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_F)+(correct_ID|plot),
                          control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                          data=peai_exp_dat)
summary(peai_perrenials)
AICc(peai_perrenials) # 644.6377

## test random intercepts
# block
peai_block<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block),
                     data=peai_exp_dat)
summary(peai_block)
AICc(peai_block) #639.632

# plot
peai_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                    data=peai_exp_dat)
summary(peai_plot)
AICc(peai_plot) #639.632

# subplot
peai_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|subplot),
                       data=peai_exp_dat)
summary(peai_subplot)
AICc(peai_subplot) #639.632

# plot within block + block
peai_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block/plot),
                          data=peai_exp_dat)
summary(peai_block_plot)
AICc(peai_block_plot) #642.1705

# subplot within plot + plot
peai_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot/subplot),
                            data=peai_exp_dat)
summary(peai_plot_subplot)
AICc(peai_plot_subplot) #642.1705

# final peai model:
peai_final<-glm.nb(developed_seeds_counted~PR_F,
                   data=peai_exp_dat)
summary(peai_final)
AICc(peai_final) #637.2174

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(peai_final, "PR_F", scale='response', rug=FALSE)

## AVBA

## test only fixed effects
avba_fixed<-glm.nb(developed_seeds_counted~PR_F,
                   data=avba_exp_dat)
summary(avba_fixed)
AICc(avba_fixed) #388.9762

## test random slopes
# number of competitors within plot
avba_comp<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                    data=avba_exp_dat)
summary(avba_comp)
AICc(avba_comp) #382.5948

# richness within plot
avba_richness<-glmer.nb(developed_seeds_counted~scale(PR_F)+(richness|plot),
                        data=avba_exp_dat)
summary(avba_richness)
AICc(avba_richness) #393.5151

# plot within block
avba_plotwblock<-glmer.nb(developed_seeds_counted~scale(PR_F)+(plot|block),
                          data=avba_exp_dat)
summary(avba_plotwblock)
AICc(avba_plotwblock) #395.8221

# presence of perrenials within plot
avba_perrenials<-glmer.nb(developed_seeds_counted~scale(PR_F)+(correct_ID|plot),
                          data=avba_exp_dat)
summary(avba_perrenials)
AICc(avba_perrenials) #396.1563 (failed to converge)

## test random intercepts
# block
avba_block<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block),
                     data=avba_exp_dat)
summary(avba_block)
AICc(avba_block) #391.3705

# plot
avba_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                    data=avba_exp_dat)
summary(avba_plot)
AICc(avba_plot) #391.332

# subplot
avba_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|subplot),
                       data=avba_exp_dat)
summary(avba_subplot)
AICc(avba_subplot) #391.3705

# plot within block + block
avba_block_plot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|block/plot),
                          data=avba_exp_dat)
summary(avba_block_plot)
AICc(avba_block_plot) #393.8431

# subplot within plot + plot
avba_plot_subplot<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot/subplot),
                            data=avba_exp_dat)
summary(avba_plot_subplot)
AICc(avba_plot_subplot) #393.8431

# final avba model:
avba_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                     data=avba_exp_dat)
summary(avba_final)
AICc(avba_final) #382.5948

# use visreg function to visualise the effect of phylogenetic relatedness on performance based on model output 
visreg(avba_final, "PR_F", scale='response', rug=FALSE)

################### FINAL MODELS

## for x axis line density plots
axis<-data.frame(x=c(0,160), y=c(0,0))

## ARCA
arca_final<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|plot),
                     data=arca_exp_dat)
summary(arca_final) #306.3924

# check model diagnostics
scatter.smooth(fitted(arca_final), residuals(arca_final, type = "pearson"),
               mgp = c(2.2, 1, 0),
               ylab = "Residuals (Pearson)",
               xlab = "Predicted on green in regulation proportion")
title("Residual plot", line = 0.7)
abline(h = 0, col="blue")

# collect confidence intervals from glmer object using emmeans and add it to the vis output for later
arca_vis_dat<-visreg(arca_final, "PR_F", scale='response', rug=FALSE)$fit
arca_em<-as.data.frame(emmeans(arca_final, specs="PR_F", type="response"))
arca_em_lower<-arca_em$response-arca_em$asymp.LCL
arca_em_upper<-arca_em$asymp.UCL-arca_em$response
for (i in 1:nrow(arca_vis_dat)){
  arca_vis_dat$visregLwr[i]<-arca_vis_dat$visregFit[i]-arca_em_lower
  arca_vis_dat$visregUpr[i]<-arca_vis_dat$visregFit[i]+arca_em_upper
}
arca_vis_dat

# use ggplot to plot visreg output (vis.out) 
arca_exp_plot <- ggplot(data = arca_vis_dat)+
  geom_line(aes(x = PR_F, y = visregFit), colour = "#CC79A7", size = 3)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=arca_exp_dat, aes(x = PR_F, y = no_inflor_with_seeds_total), size=3)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(0,28), breaks=seq(0,28,2))+
  ylab("Infloresence count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(plot.margin = unit(c(-0.05, 0.3, 0.2, 1), "cm"),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 34, colour = c("black", NA)),
        axis.text.y = element_text(size = 34, colour = c("black", NA)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length=unit(0.25, "cm"))
arca_exp_plot

# create a density plot for native range phylogenetic distance values
arca_pj_val<-arca_pj_dat$inv_ANPD
arca_native_plot <- ggdensity(arca_native_dat, "PR_F", fill="#CC79A7", alpha=0.4, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(-0.002, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = arca_pj_val), col="black", linetype = "dashed", size=3)+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2))+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.3, 0.3, -0.09, 0), "cm"),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 34, colour = "black"),
        axis.ticks.y = element_line(size = 2),
        axis.ticks.length.y=unit(0.25, "cm"))
arca_native_plot

# use ggarrange to plot experimental data and native range data together
arca_all_native<-ggarrange(arca_native_plot, arca_exp_plot, heights = c(1.5, 3.5), ncol = 1, nrow = 2)+
  theme(plot.margin = unit(c(2.5, 0.25, 0.25,3), "cm"))
arca_all_native

# MOMO
momo_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=momo_exp_dat)
summary(momo_final)

# collect confidence levels from glmer object using emmeans and add it to the vis output for later
momo_vis_dat<-visreg(momo_final, "PR_F", scale='response', rug=FALSE)$fit
momo_em<-as.data.frame(emmeans(momo_final, specs="PR_F", type="response"))
momo_em_lower<-momo_em$response-momo_em$asymp.LCL
momo_em_upper<-momo_em$asymp.UCL-momo_em$response
for (i in 1:nrow(momo_vis_dat)){
  momo_vis_dat$visregLwr[i]<-momo_vis_dat$visregFit[i]-momo_em_lower
  momo_vis_dat$visregUpr[i]<-momo_vis_dat$visregFit[i]+momo_em_upper
}
momo_vis_dat

# use ggplot to plot visreg output (vis.out) 
momo_exp_plot <- ggplot(data = momo_vis_dat)+
  geom_line(aes(x = PR_F, y = visregFit), colour = "#D55E00", size = 3)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=momo_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=3)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(0,91), breaks=seq(0,90,5))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(plot.margin = unit(c(-0.05, 0.3, 0.2, 1), "cm"),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 34, colour = c("black", NA)),
        axis.text.y = element_text(size = 34, colour = c("black", NA)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length=unit(0.25, "cm"))
momo_exp_plot

# create a density plot for native range phylogenetic distance values
momo_pj_val<-momo_pj_dat$inv_ANPD
momo_native_plot <- ggdensity(momo_native_dat, "PR_F", fill="#D55E00", alpha=0.4, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(-0.002, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = momo_pj_val), col="black", linetype = "dashed", size=3)+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"))+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.3, 0.3, -0.09, 0), "cm"),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 34, colour = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks.y = element_line(size = 2),
        axis.ticks.length.y=unit(0.25, "cm")) 
momo_native_plot 

# use ggarrange to plot experimental data and native range data together
momo_all_native<-ggarrange(momo_native_plot, momo_exp_plot, heights = c(1.5, 3.5), ncol = 1, nrow = 2)+
  theme(plot.margin = unit(c(2.5, 0.25, 0.25,3), "cm"))
momo_all_native

# PEAI
peai_final<-glm.nb(developed_seeds_counted~PR_F,
                   data=peai_exp_dat)
summary(peai_final)

# collect confidence intervals from glmer object using emmeans and add it to the vis output for later
peai_vis_dat<-visreg(peai_final, "PR_F", scale='response', rug=FALSE)$fit
peai_em<-as.data.frame(emmeans(peai_final, specs="PR_F", type="response"))
peai_em_lower<-peai_em$response-peai_em$asymp.LCL
peai_em_upper<-peai_em$asymp.UCL-peai_em$response
for (i in 1:nrow(peai_vis_dat)){
  peai_vis_dat$visregLwr[i]<-peai_vis_dat$visregFit[i]-peai_em_lower
  peai_vis_dat$visregUpr[i]<-peai_vis_dat$visregFit[i]+peai_em_upper
}
peai_vis_dat

# use ggplot to plot visreg output (vis.out) 
peai_exp_plot <- ggplot(data = peai_vis_dat)+
  geom_line(aes(x = PR_F, y = visregFit), colour = "#0072B2", size = 3)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=peai_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=3)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(0, 2400), breaks=seq(0,2400,200))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(plot.margin = unit(c(-0.05, 0.3, 0.2, 0), "cm"),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 34, colour = c("black", NA)),
        axis.text.y = element_text(size = 34, colour = c("black", NA)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length=unit(0.25, "cm"))
peai_exp_plot

# create a density plot for native range phylogenetic distance values
peai_pj_val<-peai_pj_dat$inv_ANPD
peai_native_plot <- ggdensity(peai_native_dat, "PR_F", fill="#0072B2", alpha=0.4, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(-0.002, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = peai_pj_val), col="black", linetype = "dashed", size=3)+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks.y = element_line(size = 2),
        axis.ticks.length.y=unit(0.25, "cm"))+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.3, 0.3, -0.09, 0.36), "cm"),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 34, colour = "black"))
peai_native_plot

# use ggarrange to plot experimental data and native range data together
peai_all_native<-ggarrange(peai_native_plot, peai_exp_plot, heights = c(1.5, 3.5), ncol = 1, nrow = 2)+
  theme(plot.margin = unit(c(2.5, 0.25, 0.25,3), "cm"))
peai_all_native

# AVBA
avba_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                     data=avba_exp_dat)
summary(avba_final)

# collect confidence intervals from glmer object using emmeans and add it to the vis output for later
avba_vis_dat<-visreg(avba_final, "PR_F", scale='response', rug=FALSE)$fit
avba_em<-as.data.frame(emmeans(avba_final, specs="PR_F", type="response"))
avba_em_lower<-avba_em$response-avba_em$asymp.LCL
avba_em_upper<-avba_em$asymp.UCL-avba_em$response
for (i in 1:nrow(avba_vis_dat)){
  avba_vis_dat$visregLwr[i]<-avba_vis_dat$visregFit[i]-avba_em_lower
  avba_vis_dat$visregUpr[i]<-avba_vis_dat$visregFit[i]+avba_em_upper
}
avba_vis_dat

# use ggplot to plot visreg output (vis.out) 
avba_exp_plot <- ggplot(data = avba_vis_dat)+
  geom_line(aes(x = PR_F, y = visregFit), colour = "#117733", size = 3)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=avba_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=3)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,70,5))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(plot.margin = unit(c(-0.05, 0.3, 0.2, 1.1), "cm"),
        axis.title.x = element_text(size = 40),
        axis.title.y = element_text(size = 40),
        axis.text.x = element_text(size = 34, colour = c("black", NA)),
        axis.text.y = element_text(size = 34, colour = c("black", NA)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks = element_line(size = 2),
        axis.ticks.length=unit(0.25, "cm"))
avba_exp_plot

# create a density plot for native range phylogenetic distance values
avba_pj_val<-avba_pj_dat$inv_ANPD
avba_native_plot <- ggdensity(avba_native_dat, "PR_F", fill="#117733", alpha=0.4, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,10))+
  scale_y_continuous(limits=c(-0.002, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = avba_pj_val), col="black", linetype = "dashed", size=3)+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 40),
        axis.text.y = element_text(size = 34, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect(colour = "black", size=1.2),
        axis.ticks.y = element_line(size = 2),
        axis.ticks.length.y=unit(0.25, "cm"))+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.3, 0.3, -0.09, 0.1), "cm"))
avba_native_plot

# use ggarrange to plot experimental data and native range data together
avba_all_native<-ggarrange(avba_native_plot, avba_exp_plot, heights = c(1.5, 3.5), ncol = 1, nrow = 2)+
  theme(plot.margin = unit(c(2.5, 0.25, 0.25,3), "cm"))
avba_all_native

#### plot all species together
all<-ggarrange(arca_all_native, momo_all_native,
          peai_all_native, avba_all_native,
          ncol=2, nrow=2,
          labels = c("A)", "B)", "C)", "D)"),
          font.label = list(size = 55))
all

# export/save the plot
tiff("figure_2_final.tiff", units="in", width=30, height=30, res=600)
all
dev.off()
