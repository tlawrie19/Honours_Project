# This data was used to test for problematic control plots (in association with visual assesment)

# load libraries
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggh4x)

# import data
control_testing_data<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/control_testing_data_summary.csv")
control_testing_data<-control_testing_data[which(control_testing_data$rep_metric_type=="PR_F"|
                            control_testing_data$rep_metric_type=="PR_FA"),]
control_testing_data<-control_testing_data[which(control_testing_data$plot_metric_type=="PR_F"|
                                             control_testing_data$plot_metric_type=="PR_FA"),]
control_testing_data<-na.omit(control_testing_data)

# create a column for uniq_plot
control_testing_data$foc<-NA
control_testing_data$foc[which(control_testing_data$species=="Arctotheca calendula")] <-"arca" 
control_testing_data$foc[which(control_testing_data$species=="Pentameris airoides")]<- "peai"
control_testing_data$foc[which(control_testing_data$species=="Monoculus monstrosus")]<- "momo"
control_testing_data$foc[which(control_testing_data$species=="Avena barbata")] <- "avba"

for (i in 1:nrow(control_testing_data)){
  control_testing_data$uniq_plot[i]<-paste(control_testing_data$foc[i], control_testing_data$plot[i],
                                           control_testing_data$treatment[i], sep="")
}

# create a column for plots already removed and set metric as NA
alive_list<-c("avba6thinned_control", "avba10thinned_control", 
              "peai2control", "peai2thinned_control", "peai11thinned_control")
for (i in 1:nrow(control_testing_data)){
  if (control_testing_data$uniq_plot[i] %in%alive_list==TRUE){
    control_testing_data$alive[i]<-"Y"
    control_testing_data$plot_metric_value[i]<-NA
  }
  if (control_testing_data$uniq_plot[i] %in%alive_list==FALSE){
    control_testing_data$alive[i]<-"N"
  }
}

# create a column for problematic values
prob_list<-c("arca1thinned_control", "momo10thinned_control", "avba7thinned_control")
for (i in 1:nrow(control_testing_data)){
  if (control_testing_data$uniq_plot[i] %in%prob_list==TRUE){
    control_testing_data$Outlier[i]<-"True"
  }
  if (control_testing_data$uniq_plot[i] %in%prob_list==FALSE){
    control_testing_data$Outlier[i]<-"False"
  }
}

# modify the data for plotting
control_testing_data_mod <- control_testing_data
control_testing_data_mod$plot_metric_type[which(control_testing_data_mod$plot_metric_type=="PR_F")]<-"excl. abundance"
control_testing_data_mod$plot_metric_type[which(control_testing_data_mod$plot_metric_type=="PR_FA")]<-"incl. abundance"
control_testing_data_mod$rep_metric_type[which(control_testing_data_mod$rep_metric_type=="PR_F")]<-"excl. abundance"
control_testing_data_mod$rep_metric_type[which(control_testing_data_mod$rep_metric_type=="PR_FA")]<-"incl. abundance"
control_testing_data_mod$treatment[which(control_testing_data_mod$treatment=="thinned_control")]<-"thinned control"

# plot the data
testing_plot<-na.omit(control_testing_data_mod) %>% 
  ggplot() +
  geom_boxplot(aes(x=factor(plot), y=round(rep_metric_value))) +
  geom_point(aes(x=factor(plot), y=round(plot_metric_value), colour=Outlier), pch=16, size=3.5) +
  scale_y_continuous(limits = c(0, 160), breaks=seq(0,160,40))+
  xlab("Plot") + 
  ylab("Average neighbour phylogenetic distance (myr)") +
  facet_nested(treatment+plot_metric_type~species, scale="free")+
  theme_bw()+
  theme(strip.text.x = element_text(size =28, face = "italic"),
        strip.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(size = 24),
        axis.text.y = element_text(size = 24),
        legend.title=element_text(size=28), 
        legend.text=element_text(size=28),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_color_manual(values = c("#005AB5", "#DC3220"))
testing_plot

# export the plot
tiff("figure_S2_final.tiff", units="in", width=25, height=15, res=600)
testing_plot
dev.off()
