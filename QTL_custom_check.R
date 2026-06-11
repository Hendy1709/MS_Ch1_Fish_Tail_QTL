
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

SM_data <- SM_data %>%
  select(ID,
         Standard.Length,
         Caudal.Peduncle,
         X8th.Ray,
         X3rd.Ray,
         X1st.Ray)

SM_data <- rename(SM_data,
         ID = ID,
         standard_length = Standard.Length,
         caudal_peduncle_pixels = Caudal.Peduncle,
         middle_ray_pixels = X8th.Ray,
         longest_ray_pixels = X3rd.Ray,
         outside_ray_pixels = X1st.Ray)

### Creating conversion factor
HGFRS_data <- HGFRS_data %>%
  mutate(conv.factor = standard_length/28)

SM_data <- SM_data %>%
  mutate(conv.factor = standard_length/28)

### Standardize the data
HGFRS_data_mm <- HGFRS_data %>% 
  mutate(caudal_peduncle_mm = caudal_peduncle_pixels/conv.factor,
         middle_ray_mm = middle_ray_pixels/conv.factor,
         longest_ray_mm = longest_ray_pixels/conv.factor,
         outside_ray_mm = outside_ray_pixels/conv.factor)

HGFRS_data_mm <- HGFRS_data_mm %>% 
  select(ID,
         caudal_peduncle_mm,
         middle_ray_mm,
         longest_ray_mm,
         outside_ray_mm)

SM_data_mm <- SM_data %>% 
  mutate(caudal_peduncle_mm = caudal_peduncle_pixels/conv.factor,
         middle_ray_mm = middle_ray_pixels/conv.factor,
         longest_ray_mm = longest_ray_pixels/conv.factor,
         outside_ray_mm = outside_ray_pixels/conv.factor)

SM_data_mm <- SM_data_mm %>% 
  select(ID,
         caudal_peduncle_mm,
         middle_ray_mm,
         longest_ray_mm,
         outside_ray_mm)

### implement forkedness, roundedness, and peduncle-length ratio data onto the df
#forkedness: middle ray / longest ray
#roundedness: middle ray / outside ray
#peduncle-length ratio: Caudal Peduncle / Standard Length

HGFRS_data_mm <- HGFRS_data_mm %>% 
  mutate(forkedness = longest_ray_mm/middle_ray_mm,
         roundedness = outside_ray_mm/middle_ray_mm,
         peduncle.length = caudal_peduncle_mm)

SM_data_mm <- SM_data_mm %>% 
  mutate(forkedness = longest_ray_mm/middle_ray_mm,
         roundedness = outside_ray_mm/middle_ray_mm,
         peduncle.length = caudal_peduncle_mm)

### Implement data into metadata file
HGFRS_final <- left_join(mdata, HGFRS_data_mm, by = "ID")
#Relocate columns
HGFRS_final <- HGFRS_final %>%
  relocate(caudal_peduncle_mm, middle_ray_mm, longest_ray_mm, outside_ray_mm, forkedness, roundedness, peduncle.length, .after = PC3)

SM_final <- left_join(mdata, SM_data_mm, by = "ID")
#Relocate columns
SM_final <- SM_final %>%
  relocate(caudal_peduncle_mm, middle_ray_mm, longest_ray_mm, outside_ray_mm, forkedness, roundedness, peduncle.length, .after = PC3)

###Save df as csv; for QTL mapping
write_csv(HGFRS_final, "MmAkF2.markers.HGFRS.csv", na = "")
write_csv(SM_final, "MmAkF2.markers.SM.csv", na = "")

#==========================================================
#QTL Analysis
#==========================================================

cross.mmak_HGFRS <- read.cross("csv", dir = "./",
                         file = "MmAkF2.markers.HGFRS.csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))

cross.mmak_SM <- read.cross("csv", dir = "./",
                               file = "MmAkF2.markers.SM.csv",
                               estimate.map = FALSE,
                               genotypes = c("AA","AB","BB"))

cross.mmak_HGFRS <- calc.genoprob(cross.mmak_HGFRS)
cross.mmak_SM <- calc.genoprob(cross.mmak_SM)

### Permutation testing
#Forkedness HGFRS
perm_forkedness_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness_HGFRS <- scanone(cross.mmak_HGFRS, pheno = cross.mmak_HGFRS$pheno$forkedness, method = "ehk")

#Forkedness SM
perm_forkedness_SM <- scanone(cross.mmak_SM, pheno = cross.mmak_SM$pheno$forkedness, method = "ehk", n.perm = 1000)
scanone_forkedness_SM <- scanone(cross.mmak_SM, pheno = cross.mmak_SM$pheno$forkedness, method = "ehk")

### QTL Mapping Compar
#Forkedness HGFRS
plot(scanone_forkedness_HGFRS, col = c("red"), alternate.chrid = TRUE, ylab = "LOD", ylim = c(0,3))
add.threshold(scanone_forkedness_HGFRS, perms=perm_forkedness_HGFRS, col=c("black"), lty=2)

#Forkedness SM
plot(scanone_forkedness_SM, col = c("blue"), alternate.chrid = TRUE, ylab = "LOD")
add.threshold(scanone_forkedness_SM, perms=perm_forkedness_SM, col=c("black"), lty=2)
