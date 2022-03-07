setwd("~/switchgrass/rda/tetra")
load("~/switchgrass/rda/tetra/.RData")

library(vroom)
library(vegan)
library(parallel)
library(textshape)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(tidyverse)
library(eulerr)
library(microbiome)
library(microbiomeutilities)

options(expressions = 5e5)

#### Load dataset ####

X <- vroom("rare_ASV_table_365.csv", col_names = TRUE)
ASVrows <- c("Blackwell.01.2", "Blackwell.01", "Blackwell.02", "Blackwell.03", "Blackwell.04", 
              "Blackwell.05", "Blackwell.06", "Carthage.01.2", "Carthage.01", "Carthage.02", 
              "Carthage.03", "Carthage.04", "Carthage.05", "Carthage.06", "CaveinRock.01", 
              "CaveinRock.03", "CaveinRock.04", "CaveinRock.06", "CaveinRock.07", "CaveinRock.08", 
              "Dacotah.02", "Dacotah.03", "Dacotah.04", "Dacotah.05", "Dacotah.07", "Dacotah.08", 
              "ECS.1.05", "ECS.1.06", "ECS.12.01", "ECS.12.02", "ECS.12.03", "ECS.12.05", "ECS.12.07", 
              "ECS.2.02", "ECS.2.03", "ECS.2.04", "ECS.2.06", "ECS.2.07", "ECS.10.01", "ECS.10.03", 
              "ECS.10.07", "ECS.10.08", "ECS.10.09", "ECS.10.10", "ECS.11.01", "ECS.11.02", "ECS.11.04", 
              "ECS.11.05", "ECS.11.07", "HighTide.02", "HighTide.03", "HighTide.05", "HighTide.09", 
              "HighTide.10", "Kanlow.01", "Kanlow.02", "Kanlow.03", "Kanlow.04", "Kanlow.05", "Kanlow.06", 
              "Kanlow.08", "Kanlow.10", "KY1625.01", "KY1625.02", "KY1625.03", "KY1625.04", "KY1625.06", 
              "KY1625.07", "Pathfinder.01", "Pathfinder.02", "Pathfinder.03", "Pathfinder.04", "Pathfinder.05", 
              "Pathfinder.06", "Shelter.01", "Shelter.04", "Shelter.05", "Shelter.06", "Shelter.07", "Shelter.08", 
              "Sunburst.01", "Sunburst.02", "Sunburst.04", "Sunburst.05", "Sunburst.06", "Sunburst.07", "SW102.01", 
              "SW102.02", "SW102.03", "SW102.04", "SW102.05", "SW102.06", "SW109.01", "SW109.02", "SW109.03", 
              "SW109.04", "SW109.05", "SW109.06", "SW109.07", "SW110.01", "SW110.02", "SW110.03", "SW110.04", 
              "SW110.05", "SW110.06", "SW112.01", "SW112.03", "SW112.04", "SW112.05", "SW112.06", "SW112.07", 
              "SW114.01", "SW114.02", "SW114.03", "SW114.04", "SW114.05", "SW114.07", "SW114.08", "SW115.01", 
              "SW115.02", "SW115.03", "SW115.04", "SW115.06", "SW116.02", "SW116.03", "SW116.05", "SW116.07", 
              "SW116.08", "SW116.09", "SW122.01", "SW122.02", "SW122.03", "SW122.05", "SW122.06", "SW123.01", 
              "SW123.02", "SW123.04", "SW123.05", "SW123.06", "SW123.08", "SW123.09", "SW124.01", "SW124.02", 
              "SW124.03", "SW124.04", "SW124.05", "SW124.06", "SW127.01", "SW127.03", "SW127.04", "SW127.05", 
              "SW127.06", "SW127.07", "SW128.01", "SW128.02", "SW128.03", "SW128.04", "SW128.05", "SW128.06", 
              "SW129.01", "SW129.02", "SW129.03", "SW129.04", "SW129.05", "SW129.06", "SW31.01", "SW31.03", 
              "SW31.04", "SW31.05", "SW31.06", "SW31.07", "SW33.01", "SW33.02", "SW33.03", "SW33.04", "SW33.05", 
              "SW33.06", "SW38.03", "SW38.05", "SW38.08", "SW40.01", "SW40.03", "SW40.05", "SW40.07", "SW40.08", 
              "SW40.09", "SW43.04", "SW43.05", "SW43.06", "SW43.07", "SW43.08", "SW43.09", "SW46.01", "SW46.02", 
              "SW46.04", "SW46.05", "SW46.06", "SW49.01", "SW49.03", "SW49.04", "SW49.05", "SW49.07", "SW50.01", 
              "SW50.02", "SW50.04", "SW50.05", "SW50.06", "SW51.04", "SW51.05", "SW51.06", "SW51.08", "SW51.09", 
              "SW58.01", "SW58.02", "SW58.03", "SW58.06", "SW58.07", "SW58.09", "SW63.02", "SW63.03", "SW63.04", 
              "SW63.05", "SW63.06", "SW63.07", "SW64.01", "SW64.02", "SW64.03", "SW64.04", "SW64.06", "SW64.07", 
              "SW65.03", "SW65.05", "SW65.06", "SW65.07", "SW65.08", "SW65.09", "SW65.10", "SW781.03", "SW781.05", 
              "SW781.06", "SW781.07", "SW781.08", "SW782.01", "SW782.02", "SW782.03", "SW782.09", "SW786.01", 
              "SW786.02", "SW786.03", "SW786.06", "SW786.07", "SW786.09", "SW787.01", "SW787.02", "SW787.03",
              "SW787.04", "SW787.05", "SW787.06", "SW788.01", "SW788.02", "SW788.04", "SW788.05", "SW788.06", 
              "SW788.07", "SW788.10", "SW789.01", "SW789.06", "SW789.07", "SW789.08", "SW790.01", "SW790.03", 
              "SW790.07", "SW793.02", "SW793.03", "SW793.06", "SW793.09", "SW793.10", "SW795.04", "SW795.07", 
              "SW795.09", "SW795.10", "SW796.03", "SW796.04", "SW796.05", "SW796.06", "SW796.10", "SW797.01", 
              "SW797.02", "SW797.04", "SW797.06", "SW797.09", "SW798.05", "SW798.06", "SW798.09", "SW798.10", 
              "SW799.03", "SW799.05", "SW799.07", "SW799.09", "SW799.10", "SW802.04", "SW803.03", "SW803.04", 
              "SW803.05", "SW803.06", "SW803.07", "SW803.08", "SW803.09", "SW805.01", "SW805.02", "SW805.03", 
              "SW805.04", "SW805.05", "SW805.06", "SW805.07", "SW805.08", "SW806.01", "SW806.02", "SW806.03", 
              "SW806.04", "SW806.05", "SW806.06", "SW806.08", "SW806.09", "SW808.01", "SW808.02", "SW808.03", 
              "SW808.04", "SW808.05", "SW808.06", "SW808.07", "SW809.01", "SW809.02", "SW809.03", "SW809.04", 
              "SW809.05", "SW809.06", "SW809.07", "SWG32.02", "SWG32.03", "SWG32.04", "SWG32.05", "SWG32.08", 
              "SWG39.01", "SWG39.02", "SWG39.05", "SWG39.06", "SWG39.07", "Timber.04", "Timber.06", "Timber.07", 
              "WS4U.02", "WS4U.03", "WS4U.05", "WS4U.06", "WS4U.07", "WS4U.08", "WS98SB.01", "WS98SB.02", "WS98SB.04", 
              "WS98SB.05", "WS98SB.06", "WS98SB.07", "WS98SB.09", "WS98SB.10")

X <- textshape::column_to_rownames(X, loc = 1)
row.names(X) = ASVrows


Y <- vroom("Tet_Gnum_MAF_processed.csv", col_names = TRUE) # Modified HapMap v2 SNP Matrix excluding MAF < 0.05
Y <- textshape::column_to_rownames(Y, loc = 1)


Y[1:5,1:5]
X[1:5,1:5]

Y_mat <- as.matrix(Y)
X_mat <- as.matrix(X)

X_mat <- X_mat[ rowSums(X_mat)!=0, ]
Y_mat <- subset(Y_mat, rownames(Y_mat) %in% rownames(X_mat))
X_mat <- subset(X_mat, rownames(X_mat) %in% rownames(Y_mat))

X_mat <- X_mat[ order(row.names(X_mat)), ]
Y_mat <- Y_mat[ order(row.names(Y_mat)), ]

X_mat[1:5,1:5]
Y_mat[1:5,1:5]
dim(X_mat)
dim(Y_mat)

tail(row.names(Y_mat))

#### TRANSFORM ABUNDANCE MATRIX ####
Xlog <- decostand(X_mat, "log")


#### SNP PCA ####
Y_PCA <- rda(Y_mat, scale = FALSE) #RDA with no Y matrix is just a PCA
Y_PCs <- scores(Y_PCA, display = "sites", choices=c(1:176))

summary(Y_PCA)

#### RDA #### 
sw.rda <- rda(Xlog~.,Y_PCs[,1:10], parallel = getOption("mc.cores"))

sort(sw.rda$CCA$biplot[,2])
sort(sw.rda$CCA$biplot[,1])

sw.rda

RsquareAdj(sw.rda)

# anova(sw.rda)

Z <- read.table("samdf.csv", sep = ",", header = TRUE, row.names = 1)
Z_mat <- as.matrix(Z)
Z_mat <- subset(Z_mat, rownames(Z_mat) %in% rownames(X_mat))
Z_mat[1:3,1:3]
Z_mat <- Z_mat[ order(row.names(Z_mat)), ]

dim(Z_mat)

# ggplot2, vegan plot RDA - so you can see what's happening
require(ggplot2)
require(vegan)
require(grid)
library(ggrepel)
library(plyr)
library(wesanderson)


scor = scores(sw.rda, display=c("sp", "cn", "bp"), scaling="sites", hill = F) 
scor
vif.cca(sw.rda)

barplot(sw.rda$CCA$eig/sw.rda$tot.chi, names.arg = 1:sw.rda$CCA$rank, cex.names = 0.5, ylab="Proportion of variance explained", xlab="RDA axis")

# plot in base to see what it should look like
ii = summary(sw.rda)
sp=as.data.frame(scor$species[,1:2])#Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
st=as.data.frame(ii$sites[,1:2])
yz=as.data.frame(ii$biplot[,1:2])
grp=as.data.frame(Z_mat[,"Geo"])#Grouping by Genetic Cluster
colnames(grp)="Geo"
head(grp)

#### PLOT ####

grp$Geo <- factor(grp$Geo, levels = c("North", "South", "Northeast"))
mycols <- c(West = "#C93312", South = "#096E9C", East = "#2E3436", North = "#BE3011", Northeast = "#D79C4E", ASVs = '#CDCDCD') 

ii
head(sp)
head(st)
head(yz)

sort(yz$RDA2)
yz


p <- ggplot(data = st, aes(RDA1,RDA2, color=grp$Geo, fill=grp$Geo, order=as.numeric(factor(grp$Geo), size=2))) + 
  scale_color_manual(values = c("#DCDCDC","#C93312","#D79C4E","#046C9A")) +
  geom_point() +
  stat_ellipse() +
  geom_point(data = sp, aes(RDA1,RDA2, color = "ASVs", fill = "ASVs", order = "ASVs")) +
  # geom_text_repel(data = sp, aes(RDA1, RDA2, label=row.names(sp)), size=3) + 
  theme_bw() 


p + labs(title="RDA (Showing ASVs)")


ggplot(data = sp, aes(RDA1,RDA2), size=2) +
  # geom_text(data = sp,aes(RDA1,RDA2,label=row.names(sp)),size=3) + 
  geom_point() +
  stat_ellipse(type = "norm") +
  # geom_point(data = st, aes(RDA1,RDA2, color=grp$Geo,fill=grp$Geo)) +
  theme_bw() 


#### ASV OUTLIERS ####

# Pull outlier ASVs from RDA plot.
boxplot.stats(sp$RDA1)$out
out <- boxplot.stats(sp$RDA1)$out
out_ind <- which(sp$RDA1 %in% c(out))
out_ind

boxplot.s

dim(sp)
sp_outliers1 <- sp[out_ind, ]
dim(sp_outliers1)

out <- boxplot.stats(sp$RDA2)$out
out_ind <- which(sp$RDA2 %in% c(out))
out_ind
dim(sp)
sp_outliers2 <- sp[out_ind, ]
dim(sp_outliers2)

ALL_outliers <- c(rownames(sp_outliers1), rownames(sp_outliers2))

length(ALL_outliers)

head(ALL_outliers)

# Subset Abundance Matrix by RDA Outliers
asv_outliers <- t(X_mat)
asv_outliers[1:5,1:5]

asv_outliers <- subset(asv_outliers, rownames(asv_outliers) %in% ALL_outliers) #Drops duplicates
asv_outliers <- t(asv_outliers)
dim(asv_outliers)

asv_outliers[1:5,1:5]
dim(asv_outliers)
asv_outliers <- as.data.frame(asv_outliers)

asv_outliers[1:5,1:5]
write.csv(asv_outliers, "asv_outliers.csv", row.names= TRUE)


#### phyloseq ####

library(phyloseq)
otu_table <- read.csv("asv_outliers.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)
taxa <- read.csv("TAX.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)
samdf <- read.csv("samdf.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)

taxa[1:5,1:5]
otu_table[1:5,1:5]

asvnames <- colnames(otu_table)

taxa2 <- subset(taxa, rownames(taxa) %in% asvnames)
samdf <- subset(samdf, rownames(samdf) %in% row.names(otu_table))

dim(taxa2)
dim(taxa)
dim(samdf)
dim(otu_table)
OTU <- t(data.matrix(otu_table))
TAX <- as.matrix(taxa2)

OTU[1:5,1:5]
TAX[1:5,1:5]
samdf[1:5,1:3]


ps <- phyloseq(otu_table(OTU, taxa_are_rows=TRUE), tax_table(TAX), sample_data(samdf))
ps



#### Plot Relative Abundances ####
ps@sam_data$Ecotype <- factor(ps@sam_data$Ecotype, levels = c("Lowland", "Both", "Upland"))
ps@sam_data$Ploidy <- factor(ps@sam_data$Ploidy, levels = c("Four", "Hybrid", "Eight"))
ps@sam_data$Geo <- factor(ps@sam_data$Geo, levels = c("North", "South", "Northeast"))

# colors
phylum_colors = c("#FDD5B1", "#046C9A", "#D69C4E", "#ABDDDE", "#718C85","#48706D","#7F6956","#372F21",
                  "#AC8357","#D4A05B","#BE3410","#ED713E","#A53D1C",
                  "#A45E60","#C6CDF7","#E6A0C4",
                  "#D8A499","#720E30","#C27D38","#CCC591","#798C0C",
                  "#F7E03D","#541F0F","#BD3304","#E7D2B7",
                  "#761405","#6cdd45","#655e03","#5e0fe1","#b653e9",
                  "#62c383","#7b8e7b","#a33d23","#54820a","#cfd824",
                  "#09397e","#75ceaf","#6c924d","#03de98","#ceb327",
                  "#ba8ea6")



glom <- tax_glom(ps, taxrank = "family")
glom.merge <- merge_samples(glom, "Geo")
glom.merge.prop <- transform_sample_counts(glom.merge, function(x) x / sum(x))
glom.merge.prop.melt <- psmelt(glom.merge.prop) # create dataframe from phyloseq object
# glom.merge.prop.melt.filt <- filter(glom.merge.prop.melt, Abundance > 0.02)
data_glom <- arrange(glom.merge.prop.melt, family)
data_glom

data_glom$Sample <- factor(data_glom$Sample, levels = c("North", "South", "Northeast"))
# data_glom$Sample <- factor(data_glom$Sample, levels = c("Lowland", "Both", "Upland"))
# data_glom$Sample <- factor(data_glom$Sample, levels = c("Four", "Hybrid", "Eight"))

unique(data_glom$family)

p <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill= phylum)) +
  geom_bar(stat="identity", position = "stack") +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition between Ecotypes") +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = phylum_colors) + theme(legend.position="bottom")
p


#### Venn Diagram / Core Microbiome Analysis ####
library(eulerr)
library(microbiome)
library(microbiomeutilities)

phyloGlom.f = tax_glom(ps, "family")
glomTax.f = tax_table(phyloGlom.f)[,"family"]
glomOTU.f = otu_table(phyloGlom.f, taxa_are_rows=FALSE)
glomTable.f = merge(glomOTU.f,glomTax.f,by=0,all=TRUE)
rownames(glomTable.f) = glomTable.f[,"family"]
glomTable.f$Row.names = NULL
glomTable.f$family = NULL
family.tax <- phyloGlom.f@tax_table
rownames(family.tax) <- family.tax[,"family"]
ps.family <- phyloseq(otu_table(glomTable.f, taxa_are_rows=TRUE), tax_table(family.tax), sample_data(samdf))

ps.family

table(meta(ps.family)$Geo, useNA = "always")

pseq.rel <- microbiome::transform(ps.family, "compositional")

pseq.core <- core(pseq.rel, detection = 0, prevalence = 0.1)
taxa(pseq.core)

prevs <- read.csv('snps_h2_prevs_processed.csv')
prevs$rda_outlier <- prevs$fam %in% taxa(pseq.core) 

sum(prevs$rda_outlier)/110

write_csv(prevs, 'snps_h2_prevs_rda_processed.csv')


pop <- unique(as.character(meta(pseq.rel)$Geo))
print(pop)

list_core <- c() # an empty object to store information

for (n in pop){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Geo == n) # Choose sample from Geo by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, 
                         prevalence = 0.9)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each Geo.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

mycols <- c(South = "#046C9A", North = "#FDD5B1", Northeast = "#372F21") 
plot(venn(list_core),
     fills = mycols,
     main = "Core Microbiome - Genetic Cluster", cex = 4) 


#### Tetraploid SNP PCA ####
ydf <- as.data.frame(Y_PCs)
dim(ydf)
head(grp)
dim(grp)
head(ydf)
dim(ydf)

mycols <- c(South = "#046C9A", North = "#FDD5B1", Northeast = "#D79C4E") 

theme_set(theme_light())

snp_plot <- ggplot(ydf,aes(x=PC1,y=PC2, color = grp$Geo)) + 
  geom_point() + 
  scale_color_manual(values= c("#C93312","#046C9A","#D79C4E")) +
  labs(color = "Genetic Cluster")
  
snp_plot + labs(title="Tetraploid SNP PCA")



#### VARPART ####

# pheno <- read.csv("phenotypic_measurements_2011_1.csv", row.names = 1)

p = pheno[, c('Anthesis_Date' , 'Full_Height')]
ll = read.csv("lat_lon.csv", row.names = 1)
yfam <- read.csv("family_otus.csv", row.names = 1)
yfam <- yfam[ order(row.names(yfam)), ]
tail(yfam)
dim(yfam)

yfam <- subset(yfam, rownames(yfam) %in% rownames(Z_mat))

yfam_sqrt <- sqrt(yfam)

yfam_sqrt[1:5,1:5]

dim(p)
p <- subset(p, rownames(p) %in% rownames(Z_mat))
dim(p) # phenotypes
dim(ll) #lat_lon
dim(ydf) #SNP PCs
dim(Z_mat) #Ecotype
dim(Xlog) #log-scales ASV counts
dim(X_mat) #ASV Counts

p <- p[ order(row.names(p)), ]
Z_mat <- Z_mat[ order(row.names(Z_mat)), ]
ydf <- ydf[ order(row.names(ydf)), ]
ll <- ll[ order(row.names(ll)), ]

Z_Eco <- Z_mat$Ecotype
Z_Eco <- as.data.frame(Z_Eco)
Z_Eco$Z_Eco <- as.numeric(Z_Eco$Z_Eco)



dim(Xlog)
dim(ydf)
dim(p)
dim(Z_Eco)
dim(ll)



p[is.na(p)] <- 0

ydf[,1:119]

vpS <- varpart(Y = yfam_sqrt, X = ~ ., ydf[,1:10], p, Z_Eco, data = ll) 
vpplot <- plot(vpS, digits = 4, Xnames = c('Lat_Lon','SNP PCs', 'Phenotypes', 'Ecotype'), bg = c('navy', 'tomato', 'green', 'yellow'))

total_r2 <- vpS$part$fract[length(vpS$part$fract$Adj.R.square), 3]
total_r2

#### VARPART 2 ####
tmp <- vpS$part[[2]]

vpS$part

tmp$Adj.R.square[tmp$Adj.R.square < 0] <- 0 #setting the ones less than zero to zero

tmp

vpSE <- euler(c(
  'X1'= tmp['[a] = X1 | X2+X3+X4','Adj.R.square'],
  'X2'= tmp['[b] = X2 | X1+X3+X4','Adj.R.square'],
  'X3' = tmp['[c] = X3 | X1+X2+X4','Adj.R.square'],
  'X4' = tmp['[d] = X4 | X1+X2+X3','Adj.R.square'],	
  'X1&X2' = tmp['[e]','Adj.R.square'], 
  'X3&X2' = tmp['[f]','Adj.R.square'],
  'X3&X1' = tmp['[g]','Adj.R.square'],
  'X4&X1' = tmp['[h]','Adj.R.square'],
  'X4&X2' = tmp['[i]','Adj.R.square'],
  'X4&X3' = tmp['[j]','Adj.R.square'],
  'X1&X2&X4' = tmp['[k]','Adj.R.square'],
  'X1&X2&X3' = tmp['[l]','Adj.R.square'],
  'X4&X2&X3' = tmp['[m]','Adj.R.square'],
  'X1&X4&X3' = tmp['[n]','Adj.R.square'],
  'X1&X2&X4&X3' = tmp['[o]','Adj.R.square']
), shape = 'ellipse')


sum(vpSE$fitted.values) - 0.009

vpSE

  
total_r2 <- vpS$part$fract[length(vpS$part$fract$Adj.R.square), 3]
total_r2

vpSEplot <- plot(vpSE, labels = c('Lat_Lon','SNP PCs', 'Phenotypes', 'Ecotype'))
vpSEplot



#### Family RDA ####
fams.4.rel <- microbiome::transform(ps.family, "compositional")
fams.4 <-otu_table(fams.4.rel)
fams.4 <- t(fams.4)

dim(fams.4)
dim(Y_PCs)
fams.4 <- subset(fams.4, rownames(fams.4) %in% rownames(Y_PCs))
fams.Z_mat <- subset(Z_mat, rownames(Z_mat) %in% rownames(Y_PCs))
dim(fams.4)
dim(Y_PCs)
dim(fams.Z_mat)

fams.Z_mat

fams.4 <- fams.4[,colSums(fams.4) > 0]


fams.rda <- rda(fams.4~.,Y_PCs[,1:10], parallel = getOption("mc.cores"))
fam.scor = scores(fams.rda, display=c("sp", "cn", "bp"), scaling="sites", hill = F) 

RsquareAdj(fams.rda)

fam.ii = summary(fams.rda)
fam.sp=as.data.frame(fam.scor$species[,1:2])#Depending on the drawing result, the drawing data can be enlarged or reduced to a certain extent, as follows
fam.st=as.data.frame(fam.ii$sites[,1:2])
fam.yz=as.data.frame(fam.ii$biplot[,1:2])
fam.grp=as.data.frame(fams.Z_mat[,"Geo"])#Grouping by Type
colnames(fam.grp)="Geo"
head(fam.grp)


#### PLOT ####

fam.grp$Geo <- factor(fam.grp$Geo, levels = c("North", "South", "Northeast"))


fam.ii
head(fam.sp)
head(fam.st)
head(fam.yz)

sort(fam.yz$RDA2)
fam.yz



p <- ggplot(data = fam.st, aes(RDA1,RDA2, color=fam.grp$Geo, fill=fam.grp$Geo, order=as.numeric(factor(fam.grp$Geo), size=2))) + 
  scale_color_manual(values = c("#DCDCDC","#C93312","#D79C4E","#046C9A")) +
  geom_point() +
  stat_ellipse() +
  geom_point(data = fam.sp, aes(RDA1,RDA2, color = "ASVs", fill = "ASVs", order = "ASVs")) +
  # geom_text_repel(data = sp, aes(RDA1, RDA2, label=row.names(sp)), size=3) + 
  theme_bw() 


p + labs(title="PCs/RDA relationship")


#### ASV VARPART ####

dim(X_mat) #ASV Counts
X_sqrt <- sqrt(X_mat)


dim(p) # phenotypes
dim(ll) #lat_lon
dim(ydf) #SNP PCs
dim(Z_mat) #Ecotype
dim(Xlog) #log-scales ASV counts
dim(X_mat) #ASV Counts
dim(X_sqrt)

vpS <- varpart(Y = Xlog, X = ~ ., ydf[,1:10], p, Z_Eco, data = ll) 
vpplot <- plot(vpS, digits = 4, Xnames = c('Lat_Lon','SNP PCs', 'Phenotypes', 'Ecotype'), bg = c('navy', 'tomato', 'green', 'yellow'))

total_r2 <- vpS$part$fract[length(vpS$part$fract$Adj.R.square), 3]
total_r2

#### ASV VARPART ####
tmp <- vpS$part[[2]]

vpS$part

tmp$Adj.R.square[tmp$Adj.R.square < 0] <- 0 #setting the ones less than zero to zero

tmp

vpSE <- euler(c(
  'X1'= tmp['[a] = X1 | X2+X3+X4','Adj.R.square'],
  'X2'= tmp['[b] = X2 | X1+X3+X4','Adj.R.square'],
  'X3' = tmp['[c] = X3 | X1+X2+X4','Adj.R.square'],
  'X4' = tmp['[d] = X4 | X1+X2+X3','Adj.R.square'],	
  'X1&X2' = tmp['[e]','Adj.R.square'], 
  'X3&X2' = tmp['[f]','Adj.R.square'],
  'X3&X1' = tmp['[g]','Adj.R.square'],
  'X4&X1' = tmp['[h]','Adj.R.square'],
  'X4&X2' = tmp['[i]','Adj.R.square'],
  'X4&X3' = tmp['[j]','Adj.R.square'],
  'X1&X2&X4' = tmp['[k]','Adj.R.square'],
  'X1&X2&X3' = tmp['[l]','Adj.R.square'],
  'X4&X2&X3' = tmp['[m]','Adj.R.square'],
  'X1&X4&X3' = tmp['[n]','Adj.R.square'],
  'X1&X2&X4&X3' = tmp['[o]','Adj.R.square']
), shape = 'ellipse')


sum(vpSE$fitted.values) - 0.009

vpSE


total_r2 <- vpS$part$fract[length(vpS$part$fract$Adj.R.square), 3]
total_r2

vpSEplot <- plot(vpSE, labels = c('Lat_Lon','SNP PCs', 'Phenotypes', 'Ecotype'))
vpSEplot
