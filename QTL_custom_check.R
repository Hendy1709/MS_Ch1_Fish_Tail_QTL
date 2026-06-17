
#Environment Setup===========================================================
library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(readr)
library(qtl)
library(qtlcharts)
library(dplyr)
#============================================================================

### Uploading Hendy + Flor data and Sophie's data
HGFRS_data <- read.csv("fish tail data/HGFRS.x.SM_MmAkF2_fish_tail_measurements.csv")
SM_data <- read.csv("fish tail data/sophie_data.csv")

### Upload metadata file
mdata <- read.csv("QTL_mapping_practice/Nov20_MmAk.sharedmarkers.allchrs.allind_pheno (1).csv")

### Selecting the relevant columns from each of the two main dataframes
HGFRS_data <- HGFRS_data %>%
  select(ID,
         standard_length,
         caudal_peduncle_pixels,
         middle_ray_pixels,
         longest_ray_pixels,
         outside_ray_pixels)

HGFRS_data <- rename(HGFRS_data,
                  ID = ID,
                  standard_length.HGFRS = standard_length,
                  caudal_peduncle_pixels.HGFRS = caudal_peduncle_pixels,
                  middle_ray_pixels.HGFRS = middle_ray_pixels,
                  longest_ray_pixels.HGFRS = longest_ray_pixels,
                  outside_ray_pixels.HGFRS = outside_ray_pixels)

SM_data <- SM_data %>%
  select(ID,
         Standard.Length,
         Caudal.Peduncle,
         X8th.Ray,
         X3rd.Ray,
         X1st.Ray)

SM_data <- rename(SM_data,
         ID = ID,
         standard_length.SM = Standard.Length,
         caudal_peduncle_pixels.SM = Caudal.Peduncle,
         middle_ray_pixels.SM = X8th.Ray,
         longest_ray_pixels.SM = X3rd.Ray,
         outside_ray_pixels.SM = X1st.Ray)

### Creating conversion factor
HGFRS_data <- HGFRS_data %>%
  mutate(conv.factor = standard_length.HGFRS/28)

SM_data <- SM_data %>%
  mutate(conv.factor = standard_length.SM/28)

### Standardize the data
HGFRS_data_mm <- HGFRS_data %>% 
  mutate(caudal_peduncle_mm.HGFRS = caudal_peduncle_pixels.HGFRS/conv.factor,
         middle_ray_mm.HGFRS = middle_ray_pixels.HGFRS/conv.factor,
         longest_ray_mm.HGFRS = longest_ray_pixels.HGFRS/conv.factor,
         outside_ray_mm.HGFRS = outside_ray_pixels.HGFRS/conv.factor)

HGFRS_data_mm <- HGFRS_data_mm %>% 
  select(ID,
         caudal_peduncle_mm.HGFRS,
         middle_ray_mm.HGFRS,
         longest_ray_mm.HGFRS,
         outside_ray_mm.HGFRS)

SM_data_mm <- SM_data %>% 
  mutate(caudal_peduncle_mm.SM = caudal_peduncle_pixels.SM/conv.factor,
         middle_ray_mm.SM = middle_ray_pixels.SM/conv.factor,
         longest_ray_mm.SM = longest_ray_pixels.SM/conv.factor,
         outside_ray_mm.SM = outside_ray_pixels.SM/conv.factor)

SM_data_mm <- SM_data_mm %>% 
  select(ID,
         caudal_peduncle_mm.SM,
         middle_ray_mm.SM,
         longest_ray_mm.SM,
         outside_ray_mm.SM)

### implement forkedness, roundedness, and peduncle-length ratio data onto the df
#forkedness: middle ray / longest ray
#roundedness: middle ray / outside ray
#peduncle-length ratio: Caudal Peduncle / Standard Length

HGFRS_data_mm <- HGFRS_data_mm %>% 
  mutate(forkedness.HGFRS = longest_ray_mm.HGFRS/middle_ray_mm.HGFRS,
         roundedness.HGFRS = outside_ray_mm.HGFRS/middle_ray_mm.HGFRS,
         peduncle.length.HGFRS = caudal_peduncle_mm.HGFRS)

SM_data_mm <- SM_data_mm %>% 
  mutate(forkedness.SM = longest_ray_mm.SM/middle_ray_mm.SM,
         roundedness.SM = outside_ray_mm.SM/middle_ray_mm.SM,
         peduncle.length.SM = caudal_peduncle_mm.SM)

### Implement data into metadata file
HGFRS_SM_final <- left_join(mdata, HGFRS_data_mm, by = "ID")
HGFRS_SM_final <- left_join(HGFRS_SM_final, SM_data_mm, by = "ID")

#Relocate columns
HGFRS_SM_final <- HGFRS_SM_final %>%
  relocate(caudal_peduncle_mm.HGFRS, middle_ray_mm.HGFRS, longest_ray_mm.HGFRS, outside_ray_mm.HGFRS, forkedness.HGFRS, roundedness.HGFRS, peduncle.length.HGFRS, caudal_peduncle_mm.SM, middle_ray_mm.SM, longest_ray_mm.SM, outside_ray_mm.SM, forkedness.SM, roundedness.SM, peduncle.length.SM, .after = PC3)

plot(forkedness.HGFRS ~ forkedness.SM, HGFRS_SM_final)
plot(roundedness.HGFRS ~ roundedness.SM, HGFRS_SM_final)
plot(peduncle.length.HGFRS ~ peduncle.length.SM, HGFRS_SM_final)

ids_removed <- c("409")
HGFRS_SM_final.w.removed <- HGFRS_SM_final %>%
  filter(!(ID %in% ids_removed))

plot(forkedness.HGFRS ~ forkedness.SM, HGFRS_SM_final.w.removed)
plot(roundedness.HGFRS ~ roundedness.SM, HGFRS_SM_final.w.removed)
plot(peduncle.length.HGFRS ~ peduncle.length.SM, HGFRS_SM_final.w.removed)

###Save df as csv; for QTL mapping
write_csv(HGFRS_SM_final, "MmAkF2.markers.HGFRS_and_SM.csv", na = "")


#==========================================================
#QTL Analysis
#==========================================================

cross.mmak_HGFRS_and_SM <- read.cross("csv", dir = "./",
                         file = "MmAkF2.markers.HGFRS_and_SM.csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))


cross.mmak_HGFRS_and_SM <- calc.genoprob(cross.mmak_HGFRS_and_SM)


### Permutation testing
#Forkedness
perm_forkedness_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$forkedness, method = "ehk")

perm_forkedness_SM <- scanone(cross.mmak_SM, pheno = cross.mmak_SM$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness_SM <- scanone(cross.mmak_SM, pheno = cross.mmak_SM$pheno$forkedness, method = "ehk")

#Caudal Peduncle
perm_caudped_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$peduncle.length, method = "ehk", n.perm = 1000)
scanone_caudped_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$forkedness, method = "ehk")

### QTL Mapping Compar
#Forkedness HGFRS
plot(scanone_forkedness_HGFRS, col = c("red"), alternate.chrid = TRUE, ylab = "LOD", ylim = c(0,3))
add.threshold(scanone_forkedness_HGFRS, perms=perm_forkedness_HGFRS, col=c("black"), lty=2)

#Forkedness SM
plot(scanone_forkedness_SM, col = c("blue"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_forkedness_SM, perms=perm_forkedness_SM, col=c("black"), lty=2)
