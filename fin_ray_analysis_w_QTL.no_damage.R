
#Environment Setup===========================================================
library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(readr)
library(qtl)
library(qtlcharts)
library(ggbeeswarm)
library(dplyr)
#============================================================================

#upload data
tail.fin_data.cleaned <- read.csv("fish tail data/MmAkF2_fish_tail_measurements_no_damage.csv")

#upload metadata file (includes behavior data and loci analysis)
mdata <- read.csv("QTL_mapping_practice/Nov20_MmAk.sharedmarkers.allchrs.allind_pheno (1).csv", na.strings=c("", "-"))

tail.fin_data.cleaned <- tail.fin_data.cleaned %>% #no damage
  select(ID, 
         standard_length, 
         caudal_peduncle_pixels, 
         middle_ray_pixels, 
         longest_ray_pixels,
         outside_ray_pixels)

###Creating the conversion factor
tail.fin_data.cleaned <- tail.fin_data.cleaned %>% #no damage
  mutate(conv.factor = standard_length/28)

###Standardizing the data
tail.fin_data.cleaned <- tail.fin_data.cleaned %>% #no damage
  mutate(caudal_peduncle_mm = caudal_peduncle_pixels/conv.factor,
         middle_ray_mm = middle_ray_pixels/conv.factor,
         longest_ray_mm = longest_ray_pixels/conv.factor,
         outside_ray_mm = outside_ray_pixels/conv.factor)

tail.fin_mm_cleaned <- tail.fin_data.cleaned %>% #no damage
  select(ID,
         caudal_peduncle_mm,
         middle_ray_mm,
         longest_ray_mm,
         outside_ray_mm)

### implement forkedness, roundedness, and peduncle-length ratio data onto the df
#forkedness: middle ray / longest ray
#roundedness: middle ray / outside ray
#peduncle-length ratio: Caudal Peduncle / Standard Length
tail.fin_mm_cleaned <- tail.fin_mm_cleaned %>% #no damage
  mutate(forkedness = longest_ray_mm/middle_ray_mm,
         roundedness = outside_ray_mm/middle_ray_mm,
         peduncle.length = caudal_peduncle_mm)

### Implement data into metadata file (no damage tails)
new_data.no_damage <- left_join(mdata, tail.fin_mm_cleaned, by = "ID")
#Relocate columns
new_data.no_damage <- new_data.no_damage %>%
  relocate(caudal_peduncle_mm, middle_ray_mm, longest_ray_mm, outside_ray_mm, forkedness, roundedness, peduncle.length, .after = PC3)


### Scatterplot to check
plot(middle_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage)
plot(longest_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage)
plot(outside_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage)
plot(peduncle.length ~ standard_length_mm.cleaned, data = new_data.no_damage)

### Standardizing fin ray length by standard length of fish
new_data.no_damage <- new_data.no_damage %>%
  mutate(midray_std.length = middle_ray_mm/standard_length_mm.cleaned,
         longray_std.length = longest_ray_mm/standard_length_mm.cleaned,
         outray_std.length = outside_ray_mm/standard_length_mm.cleaned,
         peduncle_std.length = peduncle.length/standard_length_mm.cleaned)

new_data.no_damage <- new_data.no_damage %>%
  relocate(midray_std.length, longray_std.length, outray_std.length, peduncle_std.length, .after = peduncle.length)

### Modelling rays against standard length to find residual models
sl.mid_cleaned <- lm(middle_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage, na.action = na.exclude) #Middle Ray residuals
res_sl.mid_cleaned <- residuals(sl.mid_cleaned)

sl.long_cleaned <- lm(longest_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage, na.action = na.exclude) #Longest Ray residuals
res_sl.long_cleaned <- residuals(sl.long_cleaned)

sl.out_cleaned <- lm(outside_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage, na.action = na.exclude) #Outside Ray residuals
res_sl.out_cleaned <- residuals(sl.out_cleaned)

sl.pedlength_cleaned <- lm(peduncle.length ~ standard_length_mm.cleaned, data = new_data.no_damage, na.action = na.exclude) #Caudal Peduncle Residuals
res_sl.pedlen_cleaned <-residuals(sl.pedlength_cleaned)

### Input residual data into df
new_data.no_damage$resid.sl_outray_cleaned <- resid(sl.out_cleaned)
new_data.no_damage$resid.sl_longray_cleaned <- resid(sl.long_cleaned)
new_data.no_damage$resid.sl_midray_cleaned <- resid(sl.mid_cleaned)
new_data.no_damage$resid.sl_pedlength_cleaned <- resid(sl.pedlength_cleaned)

new_data.no_damage <- new_data.no_damage %>% #relocate columns so that they are easier to find
  relocate(resid.sl_midray_cleaned, resid.sl_longray_cleaned, resid.sl_outray_cleaned,resid.sl_pedlength_cleaned, .after = peduncle_std.length)

### Background check on ID info
plot(middle_ray_mm ~ standard_length_mm.cleaned, data = new_data.no_damage)
lowest_ids_middle <- new_data.no_damage %>%
  slice_min(middle_ray_mm, n = 6) %>%
  pull(ID, middle_ray_mm)

lowest_ids_peduncle <- new_data.no_damage %>%
  slice_min(peduncle.length, n = 15) %>%
  select(standard_length_mm.cleaned >= 40) %>%
  pull(ID,peduncle.length)

#Save df as csv; for QTL mapping
write_csv(new_data.no_damage, "MmAkF2.all_tail_pheno.all_chr.shared_markers.no_damaged_tails.csv", na = "")

##========================================================================
## QTL Mapping
##========================================================================

cross.mmak <- read.cross("csv", dir = "./",
                         file = "MmAkF2.all_tail_pheno.all_chr.shared_markers.no_damaged_tails.csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))
cross.mmak <- calc.genoprob(cross.mmak)

#Change names; easier to track
names(cross.mmak$geno) <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10","LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", "LG22", "LG23")

### Permutation Runs
#Forkedness
perm_forkedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$forkedness, method = "ehk")

summary(perm_forkedness)
summary(scanone_forkedness)

#Roundedness
perm_roundedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$roundedness, method = "ehk", n.perm = 1000)
scanone_roundedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$roundedness, method = "ehk")

#Caudal Peduncle Depth
perm_caudal <- scanone(cross.mmak, pheno = cross.mmak$pheno$caudal_peduncle_mm, method = "ehk", n.perm = 1000)
scanone_caudal <- scanone(cross.mmak, pheno = cross.mmak$pheno$caudal_peduncle_mm, method = "ehk")

#Middle Ray Residuals
perm_midray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_midray_cleaned, method = "ehk", n.perm = 1000)
scanone_midray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_midray_cleaned, method = "ehk")

summary(perm_midray)
summary(scanone_midray)

#Longest Ray Residuals
perm_longray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_longray_cleaned, method = "ehk", n.perm = 1000)
scanone_longray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_longray_cleaned, method = "ehk")

summary(perm_longray)
summary(scanone_longray)

#Outside Ray Residuals
perm_outray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_outray_cleaned, method = "ehk", n.perm = 1000)
scanone_outray <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_outray_cleaned, method = "ehk")

summary(perm_outray)

#Caudal Peduncle Depth Residuals
perm_caudal_resid <- scanone(cross.mmak, pheno = cross.mmak$pheno$resid.sl_pedlength_cleaned, method = "ehk", n.perm = 1000)
scanone_caudal_resid <- scanone(cross.mmak, pheno = cross.mmak$pheno$caudal_peduncle_mm, method = "ehk")

summary(scanone_caudal_resid)
summary(perm_caudal_resid)

### Generate the Maps
#Middle Ray residuals
plot(scanone_midray, col = c("red"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_midray, perms=perm_midray, col=c("black"), lty=2)

#Longest Ray Residuals
plot(scanone_longray, col = c("blue"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_longray, perms=perm_longray, col=c("black"), lty=2)

#Outside Ray Residuals
plot(scanone_outray, col = c("dark green"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_outray, perms=perm_outray, col=c("black"), lty=2)

#Caudal Peduncle Depth Residuals
plot(scanone_caudal_resid, col = c("orange"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_caudal_resid, perms=perm_caudal_resid, alpha = 0.05, col=c("black"), lty=1)
add.threshold(scanone_caudal_resid, perms=perm_caudal_resid, alpha = 0.1, col=c("black"), lty=2)

#Forkedness
plot(scanone_forkedness, col = c("light blue"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_forkedness, perms=perm_forkedness, alpha = 0.05, col=c("black"), lty=1)
add.threshold(scanone_forkedness, perms=perm_forkedness, alpha = 0.1, col=c("black"), lty=2)


### Phenotypic Effect Plots=====================
new_data.no_damage <- read.csv("MmAkF2.all_tail_pheno.all_chr.shared_markers.no_damaged_tails.csv", na.strings = "")
new_data_effectplots <- new_data.no_damage[-(1:2), ]

#Caudal Peduncle Residuals LG7 
#ggplot(new_data_effectplots,aes(x = NC_036786.1_43634682, y = resid.sl_pedlength_cleaned)) +
  geom_violin(aes(fill = NC_036786.1_43634682)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle ~ Standard Length Residuals")

ggplot(new_data_effectplots,aes(x = NC_036786.1_43634682, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036786.1_43634682)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle")

ggplot(new_data_effectplots,aes(x = NC_036786.1_43634682, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036786.1_43634682)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle") +
  facet_wrap(~inferred_excluding_clear_missmatches)
  
#Caudal Peduncle Residuals LG5 
ggplot(new_data_effectplots,aes(x = NC_036784.1_2890371, y = resid.sl_pedlength_cleaned)) +
  geom_violin(aes(fill = NC_036784.1_2890371)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle ~ Standard Length Residuals")

ggplot(new_data_effectplots,aes(x = NC_036784.1_2890371, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036784.1_2890371)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle")

ggplot(new_data_effectplots,aes(x = NC_036784.1_2890371, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036784.1_2890371)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle") +
  facet_wrap(~inferred_excluding_clear_missmatches)

#Caudal Peduncle Residuals LG14
ggplot(new_data_effectplots,aes(x = NC_036793.1_8652418, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036793.1_8652418)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle")

ggplot(new_data_effectplots,aes(x = NC_036793.1_8652418, y = peduncle.length)) +
  geom_violin(aes(fill = NC_036793.1_8652418)) +
  geom_boxplot(width = 0.75, alpha = 0.75) +
  geom_beeswarm(cex = 1.2) +
  scale_x_discrete(na.translate = FALSE) +
  xlab("Genotype") +
  ylab("Caudal Peduncle") +
  facet_wrap(~inferred_excluding_clear_missmatches)
 