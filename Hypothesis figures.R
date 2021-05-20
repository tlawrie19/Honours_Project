# This code was used to create expectation of results given rejection of the null hypotheses

# load libraries
library(ggplot2)
library(ggpubr)

# import data
hyp_exp<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Moch Data/hypothesis_exp_data.csv")
hyp_exp$treatment <- factor(hyp_exp$treatment,
                       levels = c('Control','Close', 'Intermediate', 'Distant'),ordered = TRUE)
hyp_nat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Moch Data/hypothesis_nat_data.csv")

# Hypothesis 1
hyp_1<-ggplot(hyp_exp[which(hyp_exp$species=="A"),])+
  geom_boxplot(aes(x=treatment, y=fitness), fill="#C0C0C0")+
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5))+
  ylab("Fitness")+
  xlab("Treatment")+
  border() +
  theme_bw() +
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
hyp_1

# Hypothesis 2
hyp_2<-ggplot(hyp_exp[which(hyp_exp$species=="A"),])+
  geom_point(aes(x=dist, y=fitness), color="#999999")+
  geom_smooth(method="lm", aes(x=dist, y=fitness, color="#999999"), alpha=0.2)+
  scale_colour_manual(values=color_group)+
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5))+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  border() +
  theme_bw()+
  ylab("Fitness")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  rremove("legend")
hyp_2

# Hypothesis 3
color_group <- c("#999999", "#E69F00")
PJ<-ggplot(hyp_exp)+
  geom_point(aes(x=dist, y=fitness, col=species))+
  geom_smooth(method="lm", aes(x=dist, y=fitness, col=species), alpha=0.2)+
  scale_colour_manual(values=color_group)+
  scale_y_continuous(limits=c(0,36), breaks=seq(0,35,5))+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  border() +
  theme_bw()+
  ylab("Fitness")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  rremove("legend")+
  theme(plot.margin = unit(c(-0.05, 0.5, 0.5, 1.47), "cm"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
native<-ggdensity(hyp_nat, "dist", fill="species", palette=c("#999999", "#E69F00"), alpha=0.4)+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(-0.002, 0.045), breaks=seq(0, 0.04, 0.01))+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.5, 0.5, -0.09, 1.1), "cm"),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
hyp_3<-ggarrange(native, PJ, heights = c(1.5, 3.5),
                 ncol = 1, nrow = 2)+
        theme(plot.margin = unit(c(0.5, 0.5, 0.5,0.5), "cm"))
hyp_3

# create a common legend for later
leg <- get_legend(ggdensity(hyp_nat, "dist", fill="species", palette=c("#999999", "#E69F00"), alpha=0.4)+
                    geom_smooth(data=hyp_exp, method="lm", aes(x=dist, y=fitness, col=species), alpha=0)+
                    geom_point(data=hyp_exp, aes(x=dist, y=fitness, col=species))+
                    theme(legend.position="right"))

# combine all 3 plots
hyp_all<-ggarrange(ggarrange(hyp_1, hyp_2, ncol = 2, labels = c("A)", "B)"), font.label = list(size = 22)),                                              # First row with scatter plot
          ggarrange(hyp_3, legend.grob = leg, legend="right", labels = "C)", font.label = list(size = 22)),
          nrow = 2)+
          theme(plot.margin = unit(c(0, 1, 0,0), "cm"))

### export figure containing all 3 plots
tiff("figure_1_final.tiff", units="in", width=10, height=10, res=600)
hyp_all
dev.off()


############### extra for pres
# Hypothesis 3
color_group <- c("#999999", "#E69F00")
hyp_exp_trim<-hyp_exp[which(hyp_exp$species=="A"),]

PJ_orig<-ggplot(hyp_exp_trim)+
  geom_point(aes(x=dist, y=fitness), colour="#999999")+
  geom_smooth(method="lm", aes(x=dist, y=fitness), col="#999999", alpha=0.2)+
  scale_y_continuous(limits=c(0,36), breaks=seq(0,35,5))+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  border() +
  theme_bw()+
  ylab("Fitness")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  rremove("legend")+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PJ_orig

PJ_asso<-ggplot(hyp_exp_trim)+
  geom_point(aes(x=dist, y=fitness), colour="#999999")+
  geom_smooth(method="lm", aes(x=dist, y=fitness), col="#999999", alpha=0.2)+
  scale_y_continuous(limits=c(0,36), breaks=seq(0,35,5))+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  border() +
  theme_bw()+
  ylab("Fitness")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  rremove("legend")+
  theme(plot.margin = unit(c(-0.05, 0.5, 0.5, 1.47), "cm"),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
PJ_asso
hyp_nat_trim<-hyp_nat[which(hyp_nat$species=="A"),]
native_orig<-ggdensity(hyp_nat_trim, "dist", fill="#999999", alpha=0.4)+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(-0.002, 0.045), breaks=seq(0, 0.04, 0.01))+
  ylab("Density")+
  xlab("Average neighbour phylogenetic distance (myr)")+
  border() +
  theme_bw()+
  rremove("legend") + 
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
native_orig

native_asso<-ggdensity(hyp_nat_trim, "dist", fill="#999999", alpha=0.4)+
  scale_x_continuous(limits=c(0,160), breaks=seq(0,160,20))+
  scale_y_continuous(limits=c(-0.002, 0.045), breaks=seq(0, 0.04, 0.01))+
  ylab("Density")+
  border() +
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+ 
  rremove("legend") + 
  theme(plot.margin = unit(c(0.5, 0.5, -0.09, 1.1), "cm"),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 13),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
native_asso

hyp_3<-ggarrange(native_asso, PJ_asso, heights = c(1.5, 3.5),
                 ncol = 1, nrow = 2)
hyp_3

### export figure containing all 3 plots
tiff("Q31.tiff", units="in", width=6, height=6, res=600)
native_orig
dev.off()

tiff("Q32.tiff", units="in", width=6, height=6, res=600)
PJ_orig
dev.off()

tiff("Q33.tiff", units="in", width=6, height=6, res=600)
hyp_3
dev.off()
