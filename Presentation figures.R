# This code was used to create the figures used in my final seminar

# load libraries
library(ggplot2)
library(emmeans)

########### Treatment analysis
# import data
treatment_analysis<-read.csv("C:/Users/tlawr/Documents/treatment_analysis.csv")
treatment_analysis<-treatment_analysis[,1:5]

# factorise data and set levels
treatment_analysis$treatment <- factor(treatment_analysis$treatment,
                                       levels = c("Control", "Thinned control", 
                                                  "Close",
                                                  "Close-1 (exotics)", 
                                                  "Close-2 (natives)",
                                                  "Intermediate",
                                                  "Distant"),)
treatment_analysis$species <- factor(treatment_analysis$species,
                                       levels = c("arca", 
                                                  "momo",
                                                  "peai", 
                                                  "avba"))
levels(treatment_analysis$treatment) <- gsub(" ", "\n", levels(treatment_analysis$treatment))

# create data frames for each focal species
arca_dat<-treatment_analysis[which(treatment_analysis$species=="arca"),]
momo_dat<-treatment_analysis[which(treatment_analysis$species=="momo"),]
peai_dat<-treatment_analysis[which(treatment_analysis$species=="peai"),]
avba_dat<-treatment_analysis[which(treatment_analysis$species=="avba"),]

##### Plotting
### ARCA
arca_plot<-ggplot(data=arca_dat, aes(x=treatment, y=estimate, fill=Reference)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.1,
                position=position_dodge(.9))+
  xlab("Treatment") + 
  ylab("Estimate")+
  ylim(c(- 1, 1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        legend.position="top")+
  scale_fill_manual(values = c("#005AB5", "grey43"))
arca_plot
tiff("arca_treat_relat.tiff", units="in", width=8, height=7, res=600)
arca_plot
dev.off()

### momo
momo_plot<-ggplot(data=momo_dat, aes(x=treatment, y=estimate, fill=Reference)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2,
                position=position_dodge(.9))+
  xlab("Treatment") + 
  ylab("Estimate")+
  scale_y_continuous(limits = c(- 1.95, 1.3), breaks = seq(-1.5, 1, by = 0.5))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        legend.position="top")+
  scale_fill_manual(values = c("#005AB5", "grey43"))
momo_plot
tiff("momo_treat_relat.tiff", units="in", width=8, height=7, res=600)
momo_plot
dev.off()

### peai
peai_plot<-ggplot(data=peai_dat, aes(x=treatment, y=estimate, fill=Reference)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2,
                position=position_dodge(.9))+
  xlab("Treatment") + 
  ylab("Estimate")+
  scale_y_continuous(limits = c(-1.6, 1.3), breaks = seq(-1.5, 1, by = 0.5))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        legend.position="top")+
  scale_fill_manual(values = c("#005AB5", "grey43"))
peai_plot
tiff("peai_treat_relat.tiff", units="in", width=8, height=7, res=600)
peai_plot
dev.off()

### avba
avba_plot<-ggplot(data=avba_dat, aes(x=treatment, y=estimate, fill=Reference)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2,
                position=position_dodge(.9))+
  xlab("Treatment") + 
  ylab("Estimate")+
  scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.75, 0.75, by = 0.25))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        legend.position="top")+
  scale_fill_manual(values = c("#005AB5", "grey43"))
avba_plot
tiff("avba_treat_relat.tiff", units="in", width=8, height=7, res=600)
avba_plot
dev.off()

########### ANPD analysis
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

##### Plotting
## ARCA
# ANPD model
arca_final<-glmer.nb(no_inflor_with_seeds_total~scale(PR_F)+(1|plot),
                     data=arca_exp_dat)

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
  geom_line(aes(x = PR_F, y = visregFit), colour = "#CC79A7", size = 1.5)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=arca_exp_dat, aes(x = PR_F, y = no_inflor_with_seeds_total), size=1.5)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0,28), breaks=seq(0,28,2))+
  ylab("Infloresence count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())
arca_exp_plot
tiff("arca_ANPD_relat.tiff", units="in", width=8, height=7, res=600)
arca_exp_plot
dev.off()

## MOMO
# ANPD model
momo_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(1|plot),
                     control=glmerControl(optimizer="bobyqa", tolPwrss=1e-3, optCtrl=list(maxfun=2e5)),
                     data=momo_exp_dat)

# collect confidence intervals from glmer object using emmeans and add it to the vis output for later
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
  geom_line(aes(x = PR_F, y = visregFit), colour = "#D55E00", size = 1.5)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=momo_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=1.5)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0,91), breaks=seq(0,90,10))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())
momo_exp_plot
tiff("momo_ANPD_relat.tiff", units="in", width=8, height=7, res=600)
momo_exp_plot
dev.off()

## PEAI
# ANPD model
peai_final<-glm.nb(developed_seeds_counted~PR_F,
                   data=peai_exp_dat)

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
  geom_line(aes(x = PR_F, y = visregFit), colour = "#0072B2", size = 1.5)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=peai_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=1.5)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 2400), breaks=seq(0,2400,400))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())
peai_exp_plot
tiff("peai_ANPD_relat.tiff", units="in", width=8, height=7, res=600)
peai_exp_plot
dev.off()

## AVBA
# ANPD model
avba_final<-glmer.nb(developed_seeds_counted~scale(PR_F)+(no_comps_total|plot),
                     data=avba_exp_dat)

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
  geom_line(aes(x = PR_F, y = visregFit), colour = "#117733", size = 1.5)+
  geom_ribbon(aes(x= PR_F, ymin =visregLwr, ymax = visregUpr),fill = "slategrey", alpha = .175)+
  geom_point(data=avba_exp_dat, aes(x = PR_F, y = developed_seeds_counted), size=1.5)+ 
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0,70), breaks=seq(0,70,10))+
  ylab("Seed count")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw() +
  rremove("legend") + 
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())
avba_exp_plot
tiff("avba_ANPD_relat.tiff", units="in", width=8, height=7, res=600)
avba_exp_plot
dev.off()
  
#### Native range versus PJ analysis
## import data
native_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Native Range Data/native_range_plot_summary.csv")
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

####Plotting
##ARCA
arca_pj_val<-arca_pj_dat$inv_ANPD
arca_native_plot <- ggdensity(arca_native_dat, "PR_F", fill="#CC79A7", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
arca_native_plot
tiff("arca_nvi.tiff", units="in", width=8, height=6, res=600)
arca_native_plot
dev.off()

##MOMO
momo_pj_val<-momo_pj_dat$inv_ANPD
momo_native_plot <- ggdensity(momo_native_dat, "PR_F", fill="#D55E00", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = momo_pj_val), col="black", linetype = "dashed", size=1.5)+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
momo_native_plot
tiff("momo_nvi.tiff", units="in", width=8, height=6, res=600)
momo_native_plot
dev.off()

##PEAI
peai_pj_val<-peai_pj_dat$inv_ANPD
peai_native_plot <- ggdensity(peai_native_dat, "PR_F", fill="#0072B2", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  geom_vline(aes(xintercept = peai_pj_val), col="black", linetype = "dashed", size=1.5)+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
peai_native_plot
tiff("peai_nvi.tiff", units="in", width=8, height=6, res=600)
peai_native_plot
dev.off()

##AVBA
avba_pj_val<-avba_pj_dat$inv_ANPD
avba_native_plot <- ggdensity(avba_native_dat, "PR_F", fill="#117733", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.03), breaks=seq(0, 0.03, 0.005))+
  geom_vline(aes(xintercept = avba_pj_val), col="black", linetype = "dashed", size=1.5)+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
avba_native_plot
tiff("avba_nvi.tiff", units="in", width=8, height=6, res=600)
avba_native_plot
dev.off()

####Plotting without PJ line
##ARCA
arca_native_plot_npj <- ggdensity(arca_native_dat, "PR_F", fill="#CC79A7", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
arca_native_plot_npj
tiff("arca_n.tiff", units="in", width=8, height=6, res=600)
arca_native_plot_npj
dev.off()

##MOMO
momo_native_plot_npj <- ggdensity(momo_native_dat, "PR_F", fill="#D55E00", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
momo_native_plot_npj
tiff("momo_n.tiff", units="in", width=8, height=6, res=600)
momo_native_plot_npj
dev.off()

##PEAI
peai_native_plot_npj <- ggdensity(peai_native_dat, "PR_F", fill="#0072B2", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.05), breaks=seq(0, 0.05, 0.01))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
peai_native_plot_npj
tiff("peai_n.tiff", units="in", width=8, height=6, res=600)
peai_native_plot_npj
dev.off()

##AVBA
avba_native_plot_npj <- ggdensity(avba_native_dat, "PR_F", fill="#117733", alpha=0.6, col=NA)+
  coord_cartesian(xlim = c(0,160))+
  scale_x_continuous(breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(0, 0.03), breaks=seq(0, 0.03, 0.005))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(color="black", size = 11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.y = element_line(color = "black"),
        axis.line.x = element_line(color = "black"),
        panel.background = element_rect())+ 
  rremove("legend")
avba_native_plot_npj
tiff("avba_n.tiff", units="in", width=8, height=6, res=600)
avba_native_plot_npj
dev.off()
