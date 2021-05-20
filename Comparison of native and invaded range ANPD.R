# This code was used to compare ANPD values (exlcuidng abundance) between the native and invaded ranges (Q2)

# load libraries
library(ape)

# import data (angsiosperm time tree (mag) and a list of families from all the annual species in PJ (pj))
mag<-read.tree(file = "C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/Raw data/Magnoliophyta_tree.txt")
mag_list<-as.character(mag$tip.label)
pj<-read.csv(file="C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/Raw data/PJ_annual_focal_relatedness.csv")
pj_list<-as.character(unique(pj$Family))

# check matching family names
setdiff(pj_list, mag_list)

# prune angiosperm tree to only families included in PJ
angio_tree_pruned<-drop.tip(mag,mag$tip.label[-match(pj_list, mag$tip.label)])

# create a dataframe of all pariwise distances from families in the tree
pairwise_dist<-as.data.frame(round(cophenetic.phylo(angio_tree_pruned)/2))

# assign distance values for each family x focal species pair
for (i in 1:nrow(pj)){
  if ((pj$Focal[i]=="arca"|pj$Focal[i]=="momo")==TRUE){
      focal_dist<-as.data.frame(t(pairwise_dist[which(row.names(pairwise_dist)=="Asteraceae"),]))
      focal_dist$family<-as.character(rownames(focal_dist))
      colnames(focal_dist)[1]<-"distance"
      row.names(focal_dist)<-NULL
      fam<-pj$Family[i]
      fam_focal_dist<-focal_dist$distance[which(focal_dist$family==fam)]  
      pj$Dist[i]<-fam_focal_dist
  }
  if ((pj$Focal[i]=="peai"|pj$Focal[i]=="avba")==TRUE){
      focal_dist<-as.data.frame(t(pairwise_dist[which(row.names(pairwise_dist)=="Poaceae"),]))
      focal_dist$family<-as.character(rownames(focal_dist))
      colnames(focal_dist)[1]<-"distance"
      row.names(focal_dist)<-NULL
      fam<-pj$Family[i]
      fam_focal_dist<-focal_dist$distance[which(focal_dist$family==fam)]  
      pj$Dist[i]<-fam_focal_dist
  }
}

# multiple family x focal distance by number of species
pj$distxspec<-as.numeric(pj$Dist*pj$Spec)

# export PJ focal relatedness data
write.csv(pj, "C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Field Data/Raw data/PJ_annual_focal_relatedness.csv", row.names=FALSE)

# calculate total number of PJ species
pj_trimmed<-pj[which(pj$Focal=="arca"),]
pj_spec_num<-sum(pj_trimmed$Spec)
pj_spec_num

# calculate ANPD for each focal for all PJ species
focal_spec<-c("arca", "momo", "peai", "avba")
ANPD<-vector(mode="numeric", length=4)
ANPD_summ<-as.data.frame(cbind(focal_spec, ANPD))

for (i in 1:nrow(ANPD_summ)){
  foc<-ANPD_summ$focal_spec[i]
  pj_spec<-pj[which(pj$Focal==foc),]
  foc_ANPD<-round((sum(pj_spec$distxspec))/93)
  ANPD_summ$ANPD[i]<-foc_ANPD
}

# import native range data
native_dat<-read.csv("C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Native Range Data/native_range_plot_summary.csv")

## create seperate df for each species native range (with condition H + plots only <100)
arca_native_dat<-native_dat[which(native_dat$focal_species=="arca"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]
momo_native_dat<-native_dat[which(native_dat$focal_species=="momo"&native_dat$condition=="H"&native_dat$plot_size_m2<=100&native_dat$plot_size_m2>0),]
peai_native_dat<-native_dat[which(native_dat$focal_species=="peai"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]
avba_native_dat<-native_dat[which(native_dat$focal_species=="avba"&native_dat$condition=="H"&native_dat$plot_size_m2<=100),]

# create vector column for mean native range ANPD
nat_mean_ANPD<-vector(mode="numeric", length=4)
nat_mean_ANPD[1]<-round(mean(arca_native_dat$PR_F))
nat_mean_ANPD[2]<-round(mean(momo_native_dat$PR_F))
nat_mean_ANPD[3]<-round(mean(peai_native_dat$PR_F))
nat_mean_ANPD[4]<-round(mean(avba_native_dat$PR_F))
ANPD_summ<-as.data.frame(cbind(ANPD_summ, nat_mean_ANPD))
for (i in 1:nrow(ANPD_summ)){
  ANPD_summ$nat_inv_diff[i]<-ANPD_summ$nat_mean_ANPD[i]-as.numeric(ANPD_summ$ANPD[i])
}
ANPD_summ

# create vector column for minimum native range ANPD
nat_min_ANPD<-vector(mode="numeric", length=4)
nat_min_ANPD[1]<-round(min(arca_native_dat$PR_F))
nat_min_ANPD[2]<-round(min(momo_native_dat$PR_F))
nat_min_ANPD[3]<-round(min(peai_native_dat$PR_F))
nat_min_ANPD[4]<-round(min(avba_native_dat$PR_F))
ANPD_summ<-as.data.frame(cbind(ANPD_summ, nat_min_ANPD))
ANPD_summ

# create vector column for maximum native range ANPD
nat_max_ANPD<-vector(mode="numeric", length=4)
nat_max_ANPD[1]<-round(max(arca_native_dat$PR_F))
nat_max_ANPD[2]<-round(max(momo_native_dat$PR_F))
nat_max_ANPD[3]<-round(max(peai_native_dat$PR_F))
nat_max_ANPD[4]<-round(max(avba_native_dat$PR_F))
ANPD_summ<-as.data.frame(cbind(ANPD_summ, nat_max_ANPD))

# rename and change column names of the dataframe
native_v_invaded_summ<-ANPD_summ
colnames(native_v_invaded_summ)[2]<-"inv_ANPD"

# export native versus invaded ANPD data
write.csv(native_v_invaded_summ, "C:/Users/tlawr/Documents/Honours/BIOL6502 (Research)/Honours R Project/Native Range Data/ANPD_native_v_invaded_summary.csv", row.names=FALSE)
