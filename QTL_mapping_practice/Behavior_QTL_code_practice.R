library(tidyr)
library(dplyr)
library(ggplot2)
library(qtl)
library(qtlcharts)
library(ggbeeswarm)

cross.mmak <- read.cross("csv", dir = "./",
                         file = "Nov20_MmAk.sharedmarkers.allchrs.allind_pheno (1).csv",
                         estimate.map = FALSE,
                         genotypes = c("AA","AB","BB"))

cross.mmak <- calc.genoprob(cross.mmak)
names(cross.mmak$geno) <- c("LG1", "LG2", "LG3", "LG4", "LG5", "LG6", "LG7", "LG8", "LG9", "LG10","LG11", "LG12", "LG13", "LG14", "LG15", "LG16", "LG17", "LG18", "LG19", "LG20", "LG22", "LG23") 

#Permutations for straightaway speed
perm_straightaway.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$straightaway.avg.speed, method = "ehk", n.perm = 100)
scanone_straightaway.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$straightaway.avg.speed, method = "ehk")

#Permutations for center avg speed
perm_center.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$center.avg.speed, method = "ehk", n.perm = 100)
scanone_center.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$center.avg.speed, method = "ehk")

#Permutations for total avg speed
perm_total.avg.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.avg.speed, method = "ehk", n.perm = 100)
scanone_total.avg.speed <- scanone(cross.mmak, pheno = cross.mmak$pheno$total.avg.speed, method = "ehk")

summary(scanone_straightaway.speed)
summary(scanone_straightaway.speed, threshold = 3)
#summary(perm_straightaway.speed, alpha)
plot(scanone_straightaway.speed, col = c("red"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_straightaway.speed, perms=perm_straightaway.speed, col=c("black"), lty=2)
#LOD score >3 at LG 5 and LG 20

summary(scanone_center.speed)
summary(scanone_center.speed, threshold = 3)

summary(scanone_total.avg.speed)
summary(scanone_total.avg.speed, threshold = 3)

plot(scanone_straightaway.speed, scanone_center.speed , scanone_total.avg.speed , col = c("red","blue","black"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_straightaway.speed, perms=perm_straightaway.speed, col=c("black"), lty=2)

plot(scanone_straightaway.speed, scanone_center.speed , scanone_total.avg.speed ,chr = c("LG5","LG20") , col = c("red","blue","black"), alternate.chrid = TRUE, ylab = "LOD Score")
add.threshold(scanone_straightaway.speed, perms=perm_straightaway.speed, col=c("black"), lty=2)

#LG5
lodint(scanone_straightaway.speed, chr = "LG5")
bayesint(scanone_straightaway.speed, chr = "LG5")
plot(scanone_straightaway.speed, chr = "LG5", col = "red", ylab = "LOD Score")
add.threshold(scanone_straightaway.speed, perms=perm_straightaway.speed, col="black", lty=2)

#LG20
lodint(scanone_straightaway.speed, chr = "LG20")
bayesint(scanone_straightaway.speed, chr = "LG20")
plot(scanone_straightaway.speed, chr = "LG20", col = "red", ylab = "LOD Score")
add.threshold(scanone_straightaway.speed, perms=perm_straightaway.speed, col="black", lty=2)
