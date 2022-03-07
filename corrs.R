setwd("~/switchgrass/gwas_core_taxa/corrs")
load("~/switchgrass/gwas_core_taxa/corrs/.RData")

library(dplyr)
library(corrplot)
library(textshape)
library(lares)
library("ggpubr")

#### Diversity

pheno <- read.csv("phenotypic_measurements_2011_1.csv", row.names = 1)
traits <- c("Anthesis_Date", "Full_Height")
pheno <- pheno[,traits]

div <- read.csv("rare_diversity.csv", row.names = 1)

pheno <- subset(pheno, rownames(pheno) %in% rownames(fam))
div <- subset(div, rownames(div) %in% rownames(pheno))

dim(pheno)
dim(div)

div <- div[ order(row.names(div)),]
pheno <- pheno[ order(row.names(pheno)),]

data2 <- merge(div,pheno,by="row.names",all.x=TRUE)
head(data2)

data2 <- textshape::column_to_rownames(data2, loc = 1)
data2 <- as.data.frame(data2)
dim(data2)
head(data2)

samdf <- read.csv('samdf.csv')

samdf
samdf <- textshape::column_to_rownames(samdf, loc = 1)
samdf <- samdf[ order(row.names(samdf)),]
data2 <- data2[ order(row.names(data2)),]

data2 <- subset(data2, rownames(data2) %in% rownames(samdf))
samdf <- subset(samdf, rownames(samdf) %in% rownames(data2))

head(samdf)
head(data2)

data2 <- merge(samdf,data2,by="row.names",all.x=TRUE)
data2

library("ggpubr")


data3 <- subset(data2, Ecotype != "Both")
dim(data3)

data3$Geo <- factor(data3$Geo, levels = c("North", "South", "Northeast", "West", "East"))
geo_colors = c("#C93312","#096E9C","#D79C4E","#ABDDDE","#000000")
data3$Ploidy <- factor(data3$Ploidy, levels = c("Four", "Hybrid", "Eight"))

plot1 <- ggscatter(data3, x = 'Full_Height', y = 'Shannon', 
                              xlab = "Plant Height", ylab = "Shannon Diversity",
                              add = "reg.line") + scale_color_manual(values=c("#b653e9","#62c383")) +
    stat_cor(method = "pearson", label.y = 7.25)

plot2 <- ggscatter(data3, x = 'Full_Height', y = 'Shannon', color = 'Ecotype', 
          xlab = "Plant Height", ylab = "Shannon Diversity",
          add = "reg.line") + scale_color_manual(values=c("#b653e9","#62c383")) +
  facet_wrap(~Ecotype, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)

plot3 <- ggscatter(data3, x = 'Full_Height', y = 'Shannon', color = 'Ploidy', 
                                 xlab = "Plant Height", ylab = "Shannon Diversity",
                                 add = "reg.line") + scale_color_manual(values=c("#E6A0C4","#720E30","#798C0C")) +
  facet_wrap(~Ploidy, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)

plot4 <- ggscatter(data3, x = 'Full_Height', y = 'Shannon', color = 'Geo', 
          xlab = "Plant Height", ylab = "Shannon Diversity",
          add = "reg.line") + scale_color_manual(values= geo_colors) +
  facet_wrap(~Geo, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)

plot5 <- ggscatter(data3, x = 'Anthesis_Date', y = 'Shannon', 
                            xlab = "Anthesis Date", ylab = "Shannon Diversity",
                            add = "reg.line") + scale_color_manual(values=c("#b653e9", "#62c383")) +
   stat_cor(method = "pearson", label.y = 7.25)

plot6 <- ggscatter(data3, x = 'Anthesis_Date', y = 'Shannon', color = 'Ecotype', 
          xlab = "Anthesis Date", ylab = "Shannon Diversity",
          add = "reg.line") + scale_color_manual(values=c("#b653e9", "#62c383")) +
  facet_wrap(~Ecotype, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)

plot7 <- ggscatter(data3, x = 'Anthesis_Date', y = 'Shannon', color = 'Ploidy', 
                               xlab = "Anthesis Date", ylab = "Shannon Diversity",
                               add = "reg.line") + scale_color_manual(values=c("#E6A0C4","#720E30","#798C0C")) +
  facet_wrap(~Ploidy, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)

plot8 <- ggscatter(data3, x = 'Anthesis_Date', y = 'Shannon', color = 'Geo', 
          xlab = "Anthesis Date", ylab = "Shannon Diversity",
          add = "reg.line") + scale_color_manual(values= geo_colors) +
  facet_wrap(~Geo, nrow = 1) + stat_cor(method = "pearson", label.y = 7.25)


ggarrange(plot1, plot5, plot2, plot6, plot3,  plot7, plot4, plot8, labels = c("A", "D", "B", "E", "C", "F", "G", "H"), ncol = 2, nrow = 4)
