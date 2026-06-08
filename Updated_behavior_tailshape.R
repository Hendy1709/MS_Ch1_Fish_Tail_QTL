##Environment Setup===============================================================
library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(rstatix)
library(ggpubr)
##===============================================================================
##Data Wrangling=================================================================

#Upload data
sophie.data <- read.csv("fish tail data/sophie_data.csv")
#drop rows with NA
sophie.data <- drop_na(sophie.data)

#upload metadata file (includes behavior data and loci analysis)
mdata <- read.csv("Behavior_LGs4_8_17_metadata.xlsx - Metadata_andGenotypes.csv", na.strings=c("", "-"))

#Selecting the columns necessary for analysis
dat1 <- sophie.data %>%
  select(ID, 
         Standard.Length, 
         Caudal.Peduncle, 
         X1st.Ray, 
         X3rd.Ray,
         X8th.Ray)

#Changing the column names
colnames(dat1) <- c("ID",
                    "Standard.Length.pixels",
                    "Caudal.Peduncle.pixels",
                    "Outside.Ray.pixels",
                    "Longest.Ray.pixels",
                    "Middle.Ray.pixels")

#Creating the Conversion Factor
dat1 <- dat1 %>%
  mutate(conv.factor = Standard.Length.pixels/28)

#Standardizing the data
dat1 <- dat1 %>%
  mutate(Caudal.Peduncle.mm = Caudal.Peduncle.pixels/conv.factor,
         Outside.Ray.mm = Outside.Ray.pixels/conv.factor,
         Longest.Ray.mm = Longest.Ray.pixels/conv.factor,
         Middle.Ray.mm = Middle.Ray.pixels/conv.factor)

dat <- dat1 %>%
  select(ID,
         Caudal.Peduncle.mm,
         Outside.Ray.mm,
         Longest.Ray.mm,
         Middle.Ray.mm)

#implement forkedness, roundedness, and peduncle-length ratio data onto the df
#forkedness: middle ray / longest ray
#roundedness: middle ray / outside ray
#peduncle-length ratio: Caudal Peduncle / Standard Length

dat <- dat %>%
  mutate(forkedness = Middle.Ray.mm/Longest.Ray.mm,
         roundedness = Middle.Ray.mm/Outside.Ray.mm,
         std.peduncle = Caudal.Peduncle.mm)

#data wrangle the behavior data from the metadata file
behav_dat <- tibble(mdata) %>%
  select(ID, total.s.in.center, total.s.in.corner, total.avg.speed, straightaway.avg.speed)

behav_dat <- drop_na(behav_dat)
behav_dat

#combine the behavior data to dataframe
data <- left_join(dat, behav_dat, by = "ID")
data

#lots of individuals with no behavior data, should remove individuals/rows that do not have behavior data
#remove the rows with NA values in them
data <- drop_na(data)
data

#acquiring the genetic data from the metadata to merge eventually
gen_data <- mdata %>%
  select(-family,
         -P0_sire,
         -Obfamily,
         -spawn_date,
         -last_feed,
         -dissection_date,
         -time_of_dissections,
         -OB_P,
         -sex.AC,
         -P0_sire_sex,
         -mass_g,
         -standard_length_mm,
         -total.s.in.center,
         -total.s.in.corner,
         -total.avg.speed,
         -straightaway.avg.speed
  )

#merging the genetic data to the main dataframe
data <- left_join(data, gen_data, by = "ID")

##==============================================================================
##Creating most and least groups for traits=====================================

##For forkedness
# Select the top 10% of values
top_10_f <- data %>%
  slice_max(forkedness, prop = 0.1) %>%
  mutate(Group = "Most Forked")

# Select the bottom 10% of values
bottom_10_f <- data %>%
  slice_min(forkedness, prop = 0.1) %>%
  mutate(Group = "Least Forked")

#Combine the two
combined_bsa_f <- bind_rows(top_10_f, bottom_10_f)

##For Roundedness
# Select the top 10% of values
top_10_r <- data %>%
  slice_max(roundedness, prop = 0.1) %>%
  mutate(Group = "Most Rounded")

# Select the bottom 10% of values
bottom_10_r <- data %>%
  slice_min(roundedness, prop = 0.1) %>%
  mutate(Group = "Least Rounded")

#Combine the two
combined_bsa_r <- bind_rows(top_10_r, bottom_10_r)
combined_bsa_r <- combined_bsa_r %>%
  slice(-1)

##For Caudal Peduncle
##To create bulk pools
# Select the top 10% of values
top_10_p <- data %>%
  slice_max(std.peduncle, prop = 0.1) %>%
  mutate(Group = "Largest Peduncle")

# Select the bottom 10% of values
bottom_10_p <- data %>%
  slice_min(std.peduncle, prop = 0.1) %>%
  mutate(Group = "Smallest Peduncle")

#Combine the two
combined_bsa_p <- bind_rows(top_10_p, bottom_10_p)
combined_bsa_p <- combined_bsa_p %>%
  slice(-1)

##==============================================================================
##Statistical analysis==========================================================
#For Forkedness
testf1 <- wilcox.test(total.s.in.center ~ Group, data=combined_bsa_f, pair = FALSE)
testf2 <- wilcox.test(total.s.in.corner ~ Group, data=combined_bsa_f, pair = FALSE)
testf3 <- wilcox.test(total.avg.speed ~ Group, data=combined_bsa_f, pair = FALSE)
testf4 <- wilcox.test(straightaway.avg.speed ~ Group, data=combined_bsa_f, pair = FALSE)

#For Roundedness
testr1 <- wilcox.test(total.s.in.center ~ Group, data=combined_bsa_r, pair = FALSE)
testr2 <- wilcox.test(total.s.in.corner ~ Group, data=combined_bsa_r, pair = FALSE)
testr3 <- wilcox.test(total.avg.speed ~ Group, data=combined_bsa_r, pair = FALSE)
testr4 <- wilcox.test(straightaway.avg.speed ~ Group, data=combined_bsa_r, pair = FALSE)

#For Caudal Peduncle Size
testp1 <- wilcox.test(total.s.in.center ~ Group, data=combined_bsa_p, pair = FALSE)
testp2 <- wilcox.test(total.s.in.corner ~ Group, data=combined_bsa_p, pair = FALSE)
testp3 <- wilcox.test(total.avg.speed ~ Group, data=combined_bsa_p, pair = FALSE)
testp4 <- wilcox.test(straightaway.avg.speed ~ Group, data=combined_bsa_p, pair = FALSE)

##==============================================================================
##Plotting Setup================================================================

#Custom color palette
custom.color <- c(
  "Most Forked" = "#0B73E1",
  "Least Forked" = "#E09034",
  "Most Rounded" = "#E09034",
  "Least Rounded" = "#0B73E1",
  "Largest Peduncle" = "#E09034",
  "Smallest Peduncle" = "#0B73E1"
)

#Custom Theme
custom.theme <- theme_classic(base_size = 14) +  # Use larger font sizes for readability
  theme(
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.grid.major = element_blank(),
    plot.margin = margin(10, 10, 10, 10),  # Adjust plot margins for better spacing
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "none",
    plot.caption = element_text(size = 9,
                                hjust = 1,
                                vjust = 12,)
  )

##==============================================================================
##Plotting======================================================================

##Forkedness
f1 <- ggplot(combined_bsa_f, aes(x=Group, y=total.s.in.center, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.6148", y = "Time in Center") 

f2 <- ggplot(combined_bsa_f, aes(x=Group, y=total.s.in.corner, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.5643", y="Time in Corner")

f3 <- ggplot(combined_bsa_f, aes(x=Group, y=total.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.3779", y="Average Speed (mm/s)")

f4 <- ggplot(combined_bsa_f, aes(x=Group, y=straightaway.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.6565", y="Straightaway Avg. Speed (mm/s)")

#Combine plots into one for forkedness
legend_f <- get_legend(f1 + theme(legend.position = "right", legend.title=element_text(size=15), legend.text=element_text(size=10)))
fplots <- plot_grid( f1, f2, f3, f4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)  
plot_grid(fplots, legend_f, rel_widths = c(1, 0.2)) +
  draw_label("Forkedness", x=0.5, y=0, vjust=-0.5, angle= 0)


## Roundedness
r1 <- ggplot(combined_bsa_r, aes(x=Group, y=total.s.in.center, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.451", y = "Time in Center") 

r2 <- ggplot(combined_bsa_r, aes(x=Group, y=total.s.in.corner, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.8886", y="Time in Corner")

r3 <- ggplot(combined_bsa_r, aes(x=Group, y=total.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.6073", y="Average Speed (mm/s)")

r4 <- ggplot(combined_bsa_r, aes(x=Group, y=straightaway.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.9926", y="Straightaway Avg. Speed (mm/s)")

#Combine plots into one for roundedness
legend_r <- get_legend(r1 + theme(legend.position = "right", legend.title=element_text(size=15), legend.text=element_text(size=10)))
rplots <- plot_grid( r1, r2, r3, r4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)  
plot_grid(rplots, legend_r, rel_widths = c(1, 0.2)) +
  draw_label("Roundedness", x=0.5, y=0, vjust=-0.5, angle= 0)


## Caudal Peduncle
p1 <- ggplot(combined_bsa_p, aes(x=factor(Group, levels=c("Smallest Peduncle","Largest Peduncle")), y=total.s.in.center, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.9025", y = "Time in Center") 

p2 <- ggplot(combined_bsa_p, aes(x=factor(Group, levels=c("Smallest Peduncle","Largest Peduncle")), y=total.s.in.corner, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption="p: 0.7847", y="Time in Corner")

p3 <- ggplot(combined_bsa_p, aes(x=factor(Group, levels=c("Smallest Peduncle","Largest Peduncle")), y=total.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.05782", y="Average Speed (mm/s)")

p4 <- ggplot(combined_bsa_p, aes(x=factor(Group, levels=c("Smallest Peduncle","Largest Peduncle")), y=straightaway.avg.speed, fill=Group)) +
  geom_boxplot() +
  custom.theme +
  scale_fill_manual(values = custom.color) +
  labs(caption = "p: 0.4157", y="Straightaway Avg. Speed (mm/s)")

#Combine plots into one for roundedness
legend_p <- get_legend(p1 + theme(legend.position = "right", legend.title=element_text(size=15), legend.text=element_text(size=10)))
pplots <- plot_grid( p1, p2, p3, p4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)  
plot_grid(pplots, legend_p, rel_widths = c(1, 0.2)) +
  draw_label("Standardized Caudal Peduncle", x=0.5, y=0, vjust=-0.5, angle= 0)
