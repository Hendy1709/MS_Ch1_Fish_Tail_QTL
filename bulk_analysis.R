library(readr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(MASS)
library(ggbeeswarm)

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

## For forkedness===================================================================
# To create bulk pools
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


#Example using total seconds in center
lm.f1 <- lm(data = combined_bsa_f, total.s.in.center ~ forkedness)
summary(lm.f1)
anova(lm.f1)
p.f1 <- summary(lm.f1)$coefficients["forkedness", "Pr(>|t|)"]

lm.f2 <- lm(data = combined_bsa_f, total.s.in.corner ~ forkedness)
summary(lm.f2)
p.f2 <- summary(lm.f2)$coefficients["forkedness", "Pr(>|t|)"]

lm.f3 <- lm(data = combined_bsa_f, total.avg.speed ~ forkedness)
summary(lm.f3)
p.f3 <- summary(lm.f3)$coefficients["forkedness", "Pr(>|t|)"]

lm.f4 <- lm(data = combined_bsa_f, straightaway.avg.speed ~ forkedness)
summary(lm.f4)
p.f4 <- summary(lm.f4)$coefficients["forkedness", "Pr(>|t|)"]

##PLOTTING==============================================================

f1 <- ggplot(combined_bsa_f, aes(x=Group, y=total.s.in.center)) +
  geom_violin(aes(color = Group, fill = Group)) +
  ylab("Time spent in Center") +
  xlab("Forkedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.f1,3)))

ggplot(combined_bsa_f, aes(x=Group, y=total.s.in.center)) +
  geom_boxplot(aes(fill = Group)) +
  ylab("Time spent in Center") +
  xlab("Forkedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.f1,3)))

ggplot(combined_bsa_f, aes(x=Group, y=total.s.in.center, color = Group)) +
  geom_beeswarm(size = 2, cex = 2) +
  ylab("Time spent in Center") +
  xlab("Forkedness")
 


#Creating a scatterplot; forkedness against seconds in corner
f2 <- ggplot(combined_bsa_f, aes(x=forkedness, y=total.s.in.corner)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Time spent in Corners") +
  xlab("Forkedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.f2,3)))

#Creating a scatterplot; forkedness against total average speed
f3 <- ggplot(combined_bsa_f, aes(x=forkedness, y=total.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Average Speed") +
  xlab("Forkedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.f3,3)))

#Creating a scatterplot; forkedness against straightaway average speed
f4 <- ggplot(combined_bsa_f, aes(x=forkedness, y=straightaway.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Straightaway Average Speed") +
  xlab("Forkedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.f4,3)))


#Combine plots into one for forkedness
plot_grid( f1, f2, f3, f4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)

## For Caudal Peduncle===============================================================
#To create bulk pools
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

#P-values for peduncle===========
lm.p1 <- lm(data = combined_bsa_p, total.s.in.center ~ std.peduncle)
summary(lm.p1)
p.p1 <- summary(lm.p1)$coefficients["std.peduncle", "Pr(>|t|)"]

lm.p2 <- lm(data = combined_bsa_p, total.s.in.corner ~ std.peduncle)
summary(lm.p2)
p.p2 <- summary(lm.p2)$coefficients["std.peduncle", "Pr(>|t|)"]

lm.p3 <- lm(data = combined_bsa_p, total.avg.speed ~ std.peduncle)
summary(lm.p3)
p.p3 <- summary(lm.p3)$coefficients["std.peduncle", "Pr(>|t|)"]

lm.p4 <- lm(data = combined_bsa_p, straightaway.avg.speed ~ std.peduncle)
summary(lm.p4)
p.p4 <- summary(lm.p4)$coefficients["std.peduncle", "Pr(>|t|)"]

#==========Create Plots for Peduncle Length========
p1 <- ggplot(combined_bsa_p, aes(x=std.peduncle, y=total.s.in.center)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Time spent in Center") +
  xlab("Standardized Peduncle Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.p1,3)))

p2 <- ggplot(combined_bsa_p, aes(x=std.peduncle, y=total.s.in.corner)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Time spent in Corner") +
  xlab("Standardized Peduncle Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.p2,3)))

p3 <- ggplot(combined_bsa_p, aes(x=std.peduncle, y=total.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Total Average Speed") +
  xlab("Standardized Peduncle Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.p3,3)))

p4 <- ggplot(combined_bsa_p, aes(x=std.peduncle, y=straightaway.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Straightaway Average Speed") +
  xlab("Standardized Peduncle Length") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.p4,3)))

p5 <- ggplot(combined_bsa_p, aes(x=PC1, y=std.peduncle)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("PC1") +
  ylab("Standardized Peduncle") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#Combining plots into one for standardized peduncle
plot_grid(p1, p2, p3, p4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)


## For roundedness=====================================================================
#To create bulk pools
# Select the top 10% of values
top_10_r <- data %>%
  slice_max(roundedness, prop = 0.1) %>%
  mutate(Group = "More Round")

# Select the bottom 10% of values
bottom_10_r <- data %>%
  slice_min(roundedness, prop = 0.1) %>%
  mutate(Group = "Less Round")

#Combine the two
combined_bsa_r <- bind_rows(top_10_r, bottom_10_r)
combined_bsa_r <- combined_bsa_r %>%
  slice(-1)

#P-values for roundedness===========
lm.r1 <- lm(data = combined_bsa_r, total.s.in.center ~ roundedness)
summary(lm.r1)
p.r1 <- summary(lm.r1)$coefficients["roundedness", "Pr(>|t|)"]

lm.r2 <- lm(data = combined_bsa_r, total.s.in.corner ~ roundedness)
summary(lm.r2)
p.r2 <- summary(lm.r2)$coefficients["roundedness", "Pr(>|t|)"]

lm.r3 <- lm(data = combined_bsa_r, total.avg.speed ~ roundedness)
summary(lm.r3)
p.r3 <- summary(lm.r3)$coefficients["roundedness", "Pr(>|t|)"]

lm.r4 <- lm(data = combined_bsa_r, straightaway.avg.speed ~ roundedness)
summary(lm.r4)
p.r4 <- summary(lm.r4)$coefficients["roundedness", "Pr(>|t|)"]



#===Create plots for Roundedness===
#Creating a scatterplot; roundedness against seconds in center
r1 <- ggplot(combined_bsa_r, aes(x=roundedness, y=total.s.in.center)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Time spent in Center") +
  xlab("Roundedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.r1,3)))

#Creating a scatterplot; roundedness against seconds in corner
r2 <- ggplot(combined_bsa_r, aes(x=roundedness, y=total.s.in.corner)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Time spent in Corners") +
  xlab("Roundedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.r2,3)))

#Creating a scatterplot; roundedness against total average speed
r3 <- ggplot(combined_bsa_r, aes(x=roundedness, y=total.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Average Speed") +
  xlab("Roundedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.r3,3)))

#Creating a scatterplot; roundedness against straightaway average speed
r4 <- ggplot(combined_bsa_r, aes(x=roundedness, y=straightaway.avg.speed)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Straightaway Average Speed") +
  xlab("Roundedness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(caption = paste("p=",round(p.r4,3)))

#Creating a scatterplot; Principal Component Analysis for Roundedness
r5 <- ggplot(combined_bsa_r, aes(x=PC1, y=roundedness)) +
  geom_point(aes(color = Group), show.legend = FALSE, size = 2.5) +
  geom_smooth(method = "lm", color = "grey", alpha = 0.1) +
  ylab("Roundedness") +
  xlab("PC1") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#Combine plots into one for roundedness
plot_grid(r1, r2, r3, r4, labels = c("A", "B","C","D","E"), label_size=10, ncol = 2)





hist(top_10_f$forkedness)
hist(data$forkedness)
## Mann-Whitney U-test
## Use wilcox text, make sure paired = FALSE


