library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(readr)
library(qtl)
library(qtlcharts)
library(ggbeeswarm)
library(dplyr)

#Upload data
tail.fin_data <- read.csv("fish tail data/MmAkF2_fish_tail_measurements.csv")


#upload metadata file (includes behavior data and loci analysis)
mdata <- read.csv("QTL_mapping_practice/Nov20_MmAk.sharedmarkers.allchrs.allind_pheno (1).csv", na.strings=c("", "-"))

#Selecting the columns necessary for analysis
tail.fin_data <- tail.fin_data %>%
  select(ID, 
         standard_length, 
         caudal_peduncle_pixels, 
         middle_ray_pixels, 
         longest_ray_pixels,
         outside_ray_pixels)


### Creating the Conversion Factor
tail.fin_data <- tail.fin_data %>%
  mutate(conv.factor = standard_length/28)


### Standardizing the data
tail.fin_data <- tail.fin_data %>%
  mutate(caudal_peduncle_mm = caudal_peduncle_pixels/conv.factor,
         middle_ray_mm = middle_ray_pixels/conv.factor,
         longest_ray_mm = longest_ray_pixels/conv.factor,
         outside_ray_mm = outside_ray_pixels/conv.factor)


tail.fin_mm <- tail.fin_data %>%
  select(ID,
         caudal_peduncle_mm,
         middle_ray_mm,
         longest_ray_mm,
         outside_ray_mm)


### implement forkedness, roundedness, and peduncle-length ratio data onto the df
#forkedness: middle ray / longest ray
#roundedness: middle ray / outside ray
#peduncle-length ratio: Caudal Peduncle / Standard Length

tail.fin_mm <- tail.fin_mm %>%
  mutate(forkedness = longest_ray_mm/middle_ray_mm,
         roundedness = outside_ray_mm/middle_ray_mm,
         peduncle.length = caudal_peduncle_mm)


### Combine tail aspects into the metadata file
new_data <- left_join(mdata, tail.fin_mm, by = "ID")
#Relocate columns
new_data <- new_data %>%
  relocate(caudal_peduncle_mm, middle_ray_mm, longest_ray_mm, outside_ray_mm, forkedness, roundedness, peduncle.length, .after = PC3)


### Modelling and Plotting the three rays against standard length to find residuals
plot(outside_ray_mm ~ standard_length_mm.cleaned, data = new_data) #Outside ray
sl.out <- lm(outside_ray_mm ~ standard_length_mm.cleaned, data = new_data, na.action = na.exclude)
res_sl.out <- residuals(sl.out)

plot(longest_ray_mm ~ standard_length_mm.cleaned, data = new_data) #longest ray
sl.long <- lm(longest_ray_mm ~ standard_length_mm.cleaned, data = new_data, na.action = na.exclude)
res_sl.long <- residuals(sl.long)

plot(middle_ray_mm ~ standard_length_mm.cleaned, data = new_data) #middle ray
sl.mid <- lm(middle_ray_mm ~ standard_length_mm.cleaned, data = new_data, na.action = na.exclude)
res_sl.mid <- residuals(sl.mid)

### Input residual data into df
new_data$resid.sl_outray <- resid(sl.out)
new_data$resid.sl_longray <- resid(sl.long)
new_data$resid.sl_midray <- resid(sl.mid)

new_data <- new_data %>%
  relocate(resid.sl_midray, resid.sl_longray, resid.sl_outray, .after = outray_std.length)

### Standardizing fin ray length by standard length of fish
new_data <- new_data %>%
  mutate(midray_std.length = middle_ray_mm/standard_length_mm.cleaned,
         longray_std.length = longest_ray_mm/standard_length_mm.cleaned,
         outray_std.length = outside_ray_mm/standard_length_mm.cleaned)

new_data <- new_data %>%
  relocate(midray_std.length, longray_std.length, outray_std.length, .after = peduncle.length)

#Save df as csv; for QTL mapping
write_csv(new_data, "MmAkF2.all_pheno.all_chr.shared_markers.csv", na = "")

plot(forkedness ~ ID, data = new_data)
lowest_ids_fork <- new_data %>%
  slice_min(forkedness, n = 6) %>%
  pull(ID, forkedness)

plot(roundedness ~ ID, data = new_data)
highest_ids_round <- new_data %>%
  slice_max(roundedness, n = 6) %>%
  pull(ID, roundedness)

##=========================================================================
###   QTL Mapping
##=========================================================================

cross.mmak <- read.cross("csv", dir = "./",
                         file = "MmAkF2.all_pheno.all_chr.shared_markers.csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))

cross.mmak <- calc.genoprob(cross.mmak)

#Change names; easier to track
names(cross.mmak$geno) <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10","LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", "LG22", "LG23")   


### QTL for tail ratio features

#Permutations for forkedness
perm_forkedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$forkedness, method = "ehk")

#Permutations for roundedness
perm_roundedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$roundedness, method = "ehk", n.perm = 1000)
scanone_roundedness <- scanone(cross.mmak, pheno = cross.mmak$pheno$roundedness, method = "ehk")

#Permutations for caudal length
perm_caudal <- scanone(cross.mmak, pheno = cross.mmak$pheno$caudal_peduncle_mm, method = "ehk", n.perm = 1000)
scanone_caudal <- scanone(cross.mmak, pheno = cross.mmak$pheno$caudal_peduncle_mm, method = "ehk")

#QTL forkedness
plot(scanone_forkedness, col = c("red"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_forkedness, perms=perm_forkedness, col=c("black"), lty=2)

plot(scanone_forkedness, col = c("red"), alternate.chrid = TRUE, ylab = "LOD Score", chr = c("LG4")) #Look specifically at marker with significance
add.threshold(scanone_forkedness, perms=perm_forkedness, col=c("black"), lty=2)

summary(scanone_forkedness, threshold = 3)
#NC_036783.1_10094308 highest peak forkedness

#QTL roundedness
plot(scanone_roundedness, col = c("blue"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_roundedness, perms=perm_roundedness, col=c("black"), lty=2)

#QTL caudal length
plot(scanone_caudal, col = c("black"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_caudal, perms=perm_caudal, alpha=0.1, col=c("black"), lty=2)

plot(scanone_caudal, col = c("black"), alternate.chrid = TRUE, ylab = "LOD Score", chr = c("LG14"))
add.threshold(scanone_caudal, perms=perm_caudal, alpha=0.1, col=c("black"), lty=2)

summary(scanone_caudal, threshold = 3)
#NC_036793.1_7967202 highest peak caudal peduncle

### QTL for fin rays


## QTL for all 3 aspects
plot(scanone_forkedness, scanone_roundedness, scanone_caudal, col = c("red", "blue","black"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_forkedness, perms = perm_forkedness, alpha=0.1, col=c("red"), lty=2)
add.threshold(scanone_caudal, perms = perm_caudal, alpha=0.1, col=c("black"), lty=2)
add.threshold(scanone_roundedness, perms = perm_roundedness, alpha=0.1, col=c("blue"), lty=2)


### Phenotypic effect plots
new_data_effectplots <- new_data[-(1:2), ]

ggplot(subset(new_data_effectplots, !is.na(NC_036783.1_10094308)), aes(x = NC_036783.1_10094308, y = forkedness)) +
  geom_boxplot(outlier.shape = NA)

new_data_effectplots %>%
  filter(forkedness <= 1.5) %>%
  filter(forkedness >= 0.75) %>%
  filter(, !is.na(NC_036783.1_10094308)) %>%
  ggplot(aes(x = NC_036783.1_10094308, y = forkedness)) +
  geom_violin(aes(fill = NC_036783.1_10094308)) +
  geom_boxplot(width = 0.5)

new_data_effectplots %>%
  filter(, !is.na(NC_036783.1_10094308)) %>%
  ggplot(aes(x = NC_036783.1_10094308, y = forkedness)) +
  geom_violin(aes(fill = NC_036783.1_10094308)) +
  geom_boxplot()

