library(tidyr)
library(dplyr)
library(ggplot2)
library(qtl)
library(qtlcharts)
library(ggbeeswarm)

setwd("~/Dropbox/Roberts_Lab/")

cross.mmak <- read.cross("csv", dir = "./",
                         file = "Nov20_MmAk.sharedmarkers.allchrs.allind_pheno.csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))

cross.mmak <- calc.genoprob(cross.mmak)

#############length_center_instance which is total time in center############
Permutation.length_center_instance <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.s.in.center, method = "ehk", n.perm = 100)

scanone_length_center_instance <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.s.in.center, method = "ehk")

#############n_corner_instances############
Permutation.n_corner_instances <- scanone(cross.mmak, pheno = cross.mmak$pheno$instances.corner, method = "ehk", n.perm = 10000)

scanone_n_corner_instances <- scanone(cross.mmak, pheno = cross.mmak$pheno$instances.corner, method = "ehk")

#############corner_duration############
Permutation.corner_duration <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.s.in.corner, method = "ehk", n.perm = 10000)

scanone_corner_duration <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.s.in.corner, method = "ehk")


#############change in speed min1-min5############
Permutation.delta_speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$delta_speed, method = "ehk", n.perm = 10000)

scanone_delta_speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$delta_speed, method = "ehk")


#############min1_instances_slow############
Permutation.min1_slow <- scanone(cross.mmak, pheno = cross.mmak$pheno$min1.instances.slow, method = "ehk", n.perm = 10000)

scanone_min1_slow <- scanone(cross.mmak, pheno = cross.mmak$pheno$min1.instances.slow, method = "ehk")

#############straightaway_speed############
Permutation.edge_speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$straightaway.avg.speed, method = "ehk", n.perm = 100)

scanone_edge_speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$straightaway.avg.speed, method = "ehk")


#############total.distance.traveled############
Permutation.distance_traveled <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.avg.speed, method = "ehk", n.perm = 10000)

scanone_distance_traveled <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.avg.speed, method = "ehk")


#############PC1############
Permutation.PC1 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC1, method = "ehk", n.perm = 10000)

scanone_PC1 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC1, method = "ehk")

#############PC2############
Permutation.PC2 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC2, method = "ehk", n.perm = 10000)

scanone_PC2 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC2, method = "ehk")

#############PC3############
Permutation.PC3 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC3, method = "ehk", n.perm = 10000)

scanone_PC3 <- scanone(cross.mmak, pheno = cross.mmak$pheno$PC3, method = "ehk")


###########Otolith###############
Permutation.otolith <- scanone(cross.mmak, pheno = cross.mmak$pheno$Index_of_otolith_length, method = "ehk", n.perm = 10000)

scanone_otolith <- scanone(cross.mmak, pheno = cross.mmak$pheno$Index_of_otolith_length, method = "ehk")

summary(Permutation.otolith, alpha = c(.10, .05, .01))
summary(scanone_otolith)


##############sex##############
scanone_sex <- scanone(cross.mmak, pheno = cross.mmak$pheno$inferred_excluding_clear_missmatches, method = "ehk")
plot(scanone_sex)
summary(scanone_sex)

##############OB################
scanone_OB <- scanone(cross.mmak, pheno = cross.mmak$pheno$OB_num, method = "ehk")
plot(scanone_OB)

Permutation.OB <- scanone(cross.mmak, pheno = cross.mmak$pheno$OB_num, method = "ehk", n.perm = 1000)
summary(Permutation.OB, alpha = c(.10, .05, .01))
summary(scanone_OB)
plot(scanone_OB, chr = "NC_036790.1")


summary(scanone_length_center_instance)
summary(Permutation.length_center_instance, alpha = c(.10, .05, .01))
plot(scanone_length_center_instance)
#iplotScanone(scanone_length_center_instance,cross.mmak, pheno.col = 40, chr = "NC_036794.1")

summary(scanone_n_corner_instances)
summary(Permutation.n_corner_instances, alpha = c(.10, .05, .01))
#plot(scanone_n_corner_instances)

summary(scanone_corner_duration)
summary(Permutation.corner_duration, alpha = c(.10, .05, .01))
#plot(scanone_corner_duration)

summary(scanone_min1_slow)
summary(Permutation.min1_slow, alpha = c(.10, .05, .01))
#plot(scanone_min1_slow)

summary(scanone_delta_speed)
summary(Permutation.delta_speed, alpha = c(.10, .05, .01))
#plot(scanone_delta_speed)

summary(scanone_edge_speed)
summary(Permutation.edge_speed, alpha = c(.10, .05, .01))
#plot(scanone_edge_speed)

summary(scanone_distance_traveled)
summary(Permutation.distance_traveled, alpha = c(.10, .05, .01))
#plot(scanone_distance_traveled)

summary(scanone_PC1)
summary(Permutation.PC1, alpha = c(.10, .05, .01))
#plot(scanone_PC1)

summary(scanone_PC2)
summary(Permutation.PC2, alpha = c(.10, .05, .01))
#plot(scanone_PC2)

summary(scanone_PC3)
summary(Permutation.PC3, alpha = c(.10, .05, .01))
#plot(scanone_PC3)


############intervals#####################

#######LG5########
lodint(scanone_edge_speed, chr = "NC_036784.1")
bayesint(scanone_edge_speed, chr = "NC_036784.1")
plot(scanone_edge_speed, chr = "NC_036784.1") +geom_hline(yintercept = 3.93, linetype = "dashed")

lodint(scanone_n_corner_instances, chr = "NC_036784.1")
bayesint(scanone_n_corner_instances, chr = "NC_036784.1")
plot(scanone_n_corner_instances, chr = "NC_036784.1")

lodint(scanone_PC1, chr = "NC_036784.1")
bayesint(scanone_PC1, chr = "NC_036784.1")
plot(scanone_PC1, chr = "NC_036784.1")

lodint(scanone_distance_traveled, chr = "NC_036784.1")
bayesint(scanone_distance_traveled, chr = "NC_036784.1")
plot(scanone_distance_traveled, chr = "NC_036784.1")

lodint(scanone_min1_slow, chr = "NC_036784.1")
bayesint(scanone_min1_slow, chr = "NC_036784.1")
plot(scanone_min1_slow, chr = "NC_036784.1")


###########LG10##########

lodint(scanone_PC1, chr = "NC_036789.1")
bayesint(scanone_PC1, chr = "NC_036789.1")
plot(scanone_PC1, chr = "NC_036789.1")

###########LG12##########

lodint(scanone_length_center_instance, chr = "NC_036791.1")
bayesint(scanone_length_center_instance, chr = "NC_036791.1")
plot(scanone_length_center_instance, chr = "NC_036791.1")

###########LG13##########

lodint(scanone_min1_slow, chr = "NC_036792.1")
bayesint(scanone_min1_slow, chr = "NC_036792.1")
plot(scanone_min1_slow, chr = "NC_036792.1")

###########LG15##########

lodint(scanone_length_center_instance, chr = "NC_036794.1")
bayesint(scanone_length_center_instance, chr = "NC_036794.1")
plot(scanone_length_center_instance, chr = "NC_036794.1")

#######LG20########
lodint(scanone_edge_speed, chr = "NC_036799.1")
bayesint(scanone_edge_speed, chr = "NC_036799.1")
plot(scanone_edge_speed, chr = "NC_036799.1")

lodint(scanone_n_corner_instances, chr = "NC_036799.1")
bayesint(scanone_n_corner_instances, chr = "NC_036799.1")
plot(scanone_n_corner_instances, chr = "NC_036799.1")

lodint(scanone_PC1, chr = "NC_036799.1")
bayesint(scanone_PC1, chr = "NC_036799.1")
plot(scanone_PC1, chr = "NC_036799.1")

lodint(scanone_distance_traveled, chr = "NC_036799.1")
bayesint(scanone_distance_traveled, chr = "NC_036799.1")
plot(scanone_distance_traveled, chr = "NC_036799.1")

#################phenotypic effect plots#################
effect_plots <- read.csv(file = "Newmapping_MmAk_behavior_pheno_effect_plots.csv")

effect_plots$Index_of_otolith_length <- as.numeric(as.character(effect_plots$Index_of_otolith_length))
effect_plots$Average_Otolith_Length_mm_ <- as.numeric(as.character(effect_plots$Average_Otolith_Length_mm_))
effect_plots$PC1 <- as.numeric(as.character(effect_plots$PC1))
effect_plots$instances.corner <- as.numeric(as.character(effect_plots$instances.corner))
effect_plots$min1.instances.slow <- as.numeric(as.character(effect_plots$min1.instances.slow))
effect_plots$straightaway.avg.speed <- as.numeric(as.character(effect_plots$straightaway.avg.speed))
effect_plots$total.avg.speed <- as.numeric(as.character(effect_plots$total.avg.speed))
effect_plots$total.s.in.center <- as.numeric(as.character(effect_plots$total.s.in.center))

effect_plots$count_3locus_otolith <- as.character(effect_plots$count_3locus_otolith)
effect_plots$otolith_7_12 <- as.character(effect_plots$otolith_7_12)
effect_plots$otolith_7_18 <- as.character(effect_plots$otolith_7_18)
effect_plots$otolith_12_18 <- as.character(effect_plots$otolith_12_18)


Just_OB <- filter(effect_plots, OB_fam == "Y")

LG5_speed_PC1_noNAs <- filter(effect_plots, LG5_39.87733687cM_imputed != "NA")
LG5_min1slow_noNAs <- filter(effect_plots, LG5_29915931_52.17902359cM != "NA")
LG5_corner_noNAs <- filter(effect_plots, LG5_33719897_53.40448192cM_imputed != "NA")

LG7_noNAs <- filter(effect_plots, LG7_31141810_44.41366554cM_imputed_otolith != "NA")

LG10_noNAs <- filter(effect_plots, LG5_33719897_53.40448192cM_imputed != "NA")

LG12_noNAs <- filter(effect_plots, LG12_10831452_8.14452272cM_imputed != "NA")

LG13_noNAs <- filter(effect_plots, LG13_24838372_39.04458348cM_imputed != "NA")

LG18_noNAs <- filter(effect_plots, LG18_13738616_30.53935762cM_imputed != "NA")

LG15_noNAs <- filter(effect_plots, LG15_16966390_35.14312915cM_imputed != "NA")

LG20_corner_PC1_noNAs <- filter(effect_plots, LG20_26780418_52.15442856cM_imputed != "NA")
LG20_speed_noNAs <- filter(effect_plots, LG20_6382148_22.34876798cM_imputed != "NA")

otolith_3_genotypes_noNAs <- filter(effect_plots, count_3locus_otolith != "NA")


###########otolith###########
qplot(x = LG7_noNAs$LG7_31141810_44.41366554cM_imputed_otolith, y = LG7_noNAs$Index_of_otolith_length, color = LG7_noNAs$LG7_31141810_44.41366554cM_imputed_otolith) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = effect_plots$LG12_8153037_6.739341261cM_imputed, y = effect_plots$Index_of_otolith_length, color = effect_plots$LG12_8153037_6.739341261cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG18_noNAs$LG18_13738616_30.53935762cM_imputed, y = LG18_noNAs$Index_of_otolith_length, color = LG18_noNAs$LG18_13738616_30.53935762cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')

qplot(x = otolith_3_genotypes_noNAs$count_3locus_otolith, y = otolith_3_genotypes_noNAs$Index_of_otolith_length, color = otolith_3_genotypes_noNAs$count_3locus_otolith)+ geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.25) + geom_boxplot(width = 0.1) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FFEB3B", "#FDD835", "#C0CA33", "#8BC34A", "#558B2F","#0288D1", "#01579B")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = effect_plots$otolith_7_12, y = effect_plots$Index_of_otolith_length) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.25) + geom_boxplot(width = 0.1)
qplot(x = effect_plots$otolith_12_18, y = effect_plots$Index_of_otolith_length) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.25) + geom_boxplot(width = 0.1)
qplot(x = effect_plots$otolith_7_18, y = effect_plots$Index_of_otolith_length) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.25) + geom_boxplot(width = 0.1)

##########sex#############
qplot(x = effect_plots$LG7_48580168_65.44002876cM_imputed_sex, y = effect_plots$inferred_excluding_clear_missmatches) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.5) + geom_boxplot(width = 0.1)

#########OB###############
qplot(x = effect_plots$LG5_20711174_32.48231386cM_imputationNOTnecessary_Obmarker, y = effect_plots$OB_P) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.5) + geom_boxplot(width = 0.1)
qplot(x = effect_plots$LG11_25499369_15.03263121cM_imputed, y = effect_plots$OB_P) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.5) + geom_boxplot(width = 0.1)

qplot(x = Just_OB$LG5_20711174_32.48231386cM_imputationNOTnecessary_Obmarker, y = Just_OB$OB_P) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.5) + geom_boxplot(width = 0.1)

########center duration######
qplot(x = LG12_noNAs$LG12_10831452_8.14452272cM_imputed, y = LG12_noNAs$total.s.in.center, color = LG12_noNAs$LG12_10831452_8.14452272cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG15_noNAs$LG15_16966390_35.14312915cM_imputed, y = LG15_noNAs$total.s.in.center, color = LG15_noNAs$LG15_16966390_35.14312915cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')


########distance traveled######

qplot(x = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed, y = LG5_speed_PC1_noNAs$total.avg.speed, color = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG20_speed_noNAs$LG20_6382148_22.34876798cM_imputed, y = LG20_speed_noNAs$total.avg.speed, color = LG20_speed_noNAs$LG20_6382148_22.34876798cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')



########edge speed######
QTL_LG5_plot_speed <- qplot(x = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed, y = LG5_speed_PC1_noNAs$straightaway.avg.speed, color = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed) + geom_quasirandom() + geom_boxplot(width = 0.25, color = "black") + scale_color_manual(values=c(WT_color, HET_color, mutant_color)) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG20_speed_noNAs$LG20_6382148_22.34876798cM_imputed, y = LG20_speed_noNAs$straightaway.avg.speed, color = LG20_speed_noNAs$LG20_6382148_22.34876798cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')


########min1 slow######
qplot(x = LG5_min1slow_noNAs$LG5_29915931_52.17902359cM, y = LG5_min1slow_noNAs$min1.instances.slow, color = LG5_min1slow_noNAs$LG5_29915931_52.17902359cM) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG13_noNAs$LG13_24838372_39.04458348cM_imputed, y = LG13_noNAs$min1.instances.slow, color = LG13_noNAs$LG13_24838372_39.04458348cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')


########N corner instances######
qplot(x = LG5_corner_noNAs$LG5_33719897_53.40448192cM_imputed, y = LG5_corner_noNAs$instances.corner, color = LG5_corner_noNAs$LG5_33719897_53.40448192cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG20_corner_PC1_noNAs$LG20_26780418_52.15442856cM_imputed, y = LG20_corner_PC1_noNAs$instances.corner, color = LG20_corner_PC1_noNAs$LG20_26780418_52.15442856cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')


########PC1######
qplot(x = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed, y = LG5_speed_PC1_noNAs$PC1, color = LG5_speed_PC1_noNAs$LG5_39.87733687cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG10_noNAs$LG5_33719897_53.40448192cM_imputed, y = LG10_noNAs$PC1, color = LG10_noNAs$LG5_33719897_53.40448192cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')
qplot(x = LG20_corner_PC1_noNAs$LG20_26780418_52.15442856cM_imputed, y = LG20_corner_PC1_noNAs$PC1, color = LG20_corner_PC1_noNAs$LG20_26780418_52.15442856cM_imputed) + geom_violin(alpha = 0.5) + geom_quasirandom(width = 0.4) + geom_boxplot(width = 0.25) + scale_color_manual(values=c("#FDD835","#4CAF50","#2196F3")) + theme_minimal() + labs(color = 'Genotype')






