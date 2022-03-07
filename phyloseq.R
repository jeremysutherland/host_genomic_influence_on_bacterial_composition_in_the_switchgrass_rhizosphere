setwd("~/switchgrass/16s/dada2/csvs")
load("~/switchgrass/16s/dada2/csvs/.RData")

#### Phyloseq ####

library(phyloseq)
library(vegan)
library(Biostrings)
library(ggplot2)
library(wesanderson)
library(car)
library(tidyr)
library(dplyr)
library(plyr)
library(vroom)
theme_set(theme_light())

otu_table <- read.csv("otu_table_reheader.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1)
taxa <- read.csv("taxa.csv", header=T, dec=".",sep=",", skipNul = T)
samdf <- read.csv("merge.csv", header=T, dec=".",sep=",", skipNul = T, row.names = 1) # meta data table generated from https://doi.org/10.1371/journal.pgen.1003215 Supp. Table 1

rownames(taxa) <- row.names(otu_table) 
OTU <- data.matrix(otu_table)
TAX <- as.matrix(taxa)

OTU[1:5,1:5]
TAX[1:5,1:5]
samdf[1:5,1:5]


ps <- phyloseq(otu_table(OTU, taxa_are_rows=TRUE), tax_table(TAX), sample_data(samdf))
ps

#### Prune ####
#Subset taxa
nrow(ps@tax_table)
ps.prune = subset_taxa(ps, domain=="Bacteria")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, phylum !="NA")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, family != "Mitochondria")
nrow(ps.prune@tax_table)
ps.prune = subset_taxa(ps.prune, order != "Chloroplast")
nrow(ps.prune@tax_table)

# prune OTUs that are not present in at least one sample
nrow(ps.prune@otu_table)
ps.prune <- prune_taxa(taxa_sums(ps.prune) > 0, ps.prune)
nrow(ps.prune@otu_table)

# prune samples with no OTUS
ncol(ps.prune@otu_table)
ps.prune <- prune_samples(sample_sums(ps.prune) > 0, ps.prune)
ncol(ps.prune@otu_table)

#### Rarefy ####
min(colSums(ps.prune@otu_table))
max(colSums(ps.prune@otu_table))

# rarecurve(t(otu_table(ps.prune)), sample = ,step=50, cex=0.5, label = F) 
ps.rare <- rarefy_even_depth(ps.prune, sample.size = 1000, rngseed = 336, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
ps.rare <- prune_taxa(taxa_sums(ps.rare) > 0, ps.rare)
ps.rare <- prune_samples(sample_sums(ps.rare) > 0, ps.rare)

min(colSums(ps.rare@otu_table))
max(colSums(ps.rare@otu_table))

ps
ps.prune
ps.rare

domain = get_taxa_unique(ps.rare, "domain")
phylum = get_taxa_unique(ps.rare, "phylum")
classes = get_taxa_unique(ps.rare, "class")
orders = get_taxa_unique(ps.rare, "order")
family = get_taxa_unique(ps.rare, "family")
genus = get_taxa_unique(ps.rare, "genus")
species = get_taxa_unique(ps.rare, "species")

#These should = zero
sum(genus == "Mitochondria", na.rm = TRUE)
sum(family == "Mitochondria", na.rm = TRUE)
sum(orders == "Mitochondria", na.rm = TRUE)
sum(classes == "Mitochondria", na.rm = TRUE)
sum(phylum == "Mitochondria", na.rm = TRUE)
sum(domain == "Mitochondria", na.rm = TRUE)
sum(genus == "Chloroplast", na.rm = TRUE)
sum(family == "Chloroplast", na.rm = TRUE)
sum(orders == "Chloroplast", na.rm = TRUE)
sum(classes == "Chloroplast", na.rm = TRUE)
sum(phylum == "Chloroplast", na.rm = TRUE)
sum(domain == "Chloroplast", na.rm = TRUE)



####Transform####
# Relative Abundance Transformation
ps.even = transform_sample_counts(ps.rare, function(otu) otu/sum(otu))

ps.even
ps.rare


length(get_taxa_unique(ps.rare, "domain"))
length(get_taxa_unique(ps.rare, "phylum")) 
length(get_taxa_unique(ps.rare, "class"))
length(get_taxa_unique(ps.rare, "order"))
length(get_taxa_unique(ps.rare, "family")) 
length(get_taxa_unique(ps.rare, "genus"))
length(get_taxa_unique(ps.rare, "species"))

#### Ordinations ####

### PCoA

ps.even@sam_data$Ecotype <- factor(ps.even@sam_data$Ecotype, levels = c("Lowland", "Both", "Upland"))
ps.even@sam_data$Ploidy <- factor(ps.even@sam_data$Ploidy, levels = c("Four", "Hybrid", "Eight"))
ps.even@sam_data$Geo <- factor(ps.even@sam_data$Geo, levels = c("North", "South", "Northeast", "West", "East"))

#All Genotypes
ps.pcoa <- ordinate(ps.even, "PCoA", "bray")
ps.pcoa

#Tetraploids and Octoploids seperately
ps.even.4 <- subset_samples(ps.even, Ploidy == "Four")
ps.even.8 <- subset_samples(ps.even, Ploidy == "Eight")
ps.pcoa.4 <- ordinate(ps.even.4, "PCoA", "bray")
ps.pcoa.8 <- ordinate(ps.even.8, "PCoA", "bray")

write.csv(ps.pcoa$vectors, "vectors.csv")
write.csv(ps.pcoa.4$vectors, "4vectors.csv")
write.csv(ps.pcoa.8$vectors, "8vectors.csv")

eco.pcoa <- plot_ordination(ps.even, ps.pcoa, type="sample", color= "Ecotype", title="PCoA - Ecotype") 
eco.pcoa + geom_point(size = 2) + stat_ellipse(type = "norm", linetype = 1) + scale_color_manual(values = c("#D69C4E","#000000","#046C9A")) + theme_bw()

geo.pcoa <- plot_ordination(ps.even.4, ps.pcoa.4, type="sample", color= "Geo")#, title="PCoA - Genetic Cluster") 
geo.pcoa + geom_point(size = 2) + scale_color_manual(values = c("#C93312","#046C9A","#000000"), name = "Genetic Cluster") + stat_ellipse(type = "norm", linetype = 1) + theme_bw()

ploidy.pcoa <- plot_ordination(ps.even, ps.pcoa, type="sample", color= "Ploidy", title="PCoA - Ploidy")
ploidy.pcoa + geom_point(size = 2) + stat_ellipse(type = "norm", linetype = 1) + scale_color_manual(values = c("#D69C4E","#000000","#046C9A")) + theme_bw()

pop.pcoa <- plot_ordination(ps.even, ps.pcoa, type="sample", color= "Population", title="PCoA - Population")
pop.pcoa #+ geom_point(size = 2) + stat_ellipse(type = "norm", linetype = 1) + scale_color_manual(values = rev(wes_palette("Darjeeling2", n=5))) + theme_bw()


#### Boxplots ####
library(ggpubr)
library(cowplot)
# # http://statweb.stanford.edu/~susan/summer12/phyloseq-demosh.html
ps.rare@sam_data$Ecotype <- factor(ps.rare@sam_data$Ecotype, levels = c("Lowland", "Both", "Upland"))
ps.rare@sam_data$Ploidy <- factor(ps.rare@sam_data$Ploidy, levels = c("Four", "Hybrid", "Eight"))
ps.rare@sam_data$Geo <- factor(ps.rare@sam_data$Geo, levels = c("North", "South", "Northeast", "West", "East"))

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


p1 <- plot_richness(ps.rare, x="Ecotype", measures=c("Chao1", "Shannon"), color = "Ecotype") + 
  geom_violin() + 
  geom_point (size = 4, alpha = 0.3, position = position_jitter(width = 0.3)) + 
  scale_color_manual(values = c(phylum_colors[30],phylum_colors[7],phylum_colors[31])) 

p1 +
  theme(axis.text=element_text(size=16, face = "bold"), axis.title=element_text(size=16)) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.title.y = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(plot.title = element_text(size=22))

p2 <- plot_richness(ps.rare, x="Ploidy", measures=c("Chao1", "Shannon"), color = "Ploidy") + 
  geom_violin() + 
  geom_point (size = 3, alpha = 0.3, position = position_jitter(width = 0.3)) + 
  scale_color_manual(values = c(phylum_colors[16],phylum_colors[18],phylum_colors[21])) 

p2 +
  theme(axis.text=element_text(size=16, face = "bold"), axis.title=element_text(size=16)) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.title.y = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(plot.title = element_text(size=22))

p3 <- plot_richness(ps.rare, x="Geo", measures=c("Chao1", "Shannon"), color = "Geo") + 
  geom_violin() + geom_point (size = 3, alpha = 0.3, position = position_jitter(width = 0.3)) + 
  scale_color_manual(values = c("#C93312","#046C9A","#D69C4E","#ABDDDE","#000000"), name = "Genetic Cluster")

p3 +
  labs(x = "Genetic Cluster") +
  theme(axis.text=element_text(size=16, face = "bold"), axis.title=element_text(size=16)) +
  theme(axis.title.x = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.title.y = element_text(margin = margin(t = 20, r = 20, b = 20, l = 20))) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  theme(plot.title = element_text(size=22))



ggdraw() +
  draw_plot(p1, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(p2, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(p3, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

#### Statistics ####
ps.rare
div <- estimate_richness(ps.rare, split = TRUE)
dim(div)
write.csv(div, "diversity.csv")
shannon <- estimate_richness(ps.rare, split = TRUE, measures = "Shannon")
shannon

df_S <- data.frame(sample=names(shannon), Diversity = shannon, measure =rep("Shannon", length(shannon)))
metadat_d <- cbind(df_S, ps.even@sam_data)

metadat_d

leveneTest(Shannon ~ Ploidy , data=metadat_d)

#Permanova test using adonis function
dist = sqrt(phyloseq::distance(ps.even, "bray"))
distm = as.matrix(dist)

heatmap(distm, scale="column")

dist = sqrt(phyloseq::distance(ps.even, "bray"))
distnames <- names(dist)
distnames

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(ps.even))
sampledf

samdf2 <- subset(samdf, rownames(samdf) %in% distnames)
samdf2 <- samdf2[match(distnames, rownames(samdf2)),]
samdf2

# NPMANOVA test all groups
adonis2(dist ~ Ecotype * Ploidy * Geo, data = samdf2, method = "euclidean",
       permutations = 999)


# NPMANOVA test Ecotypes
adonis2(dist ~ Ecotype, data = samdf2, method = "euclidean",
       permutations = 999)

# NPMANOVA test Ploidy
adonis(dist ~ Ploidy, data = sampledf, method = "euclidean",
       permutations = 999) 

# NPMANOVA test Genetic Clusters
adonis(dist ~ Geo, data = sampledf, method = "euclidean",
       permutations = 999) 


rda(dist) # rda w/out constraining matrix = pcoa

samdf2 <- subset(samdf, rownames(samdf) %in% distnames)
samdf2 <- samdf2[match(distnames, rownames(samdf2)),]
samdf2

nrow(samdf2)
nrow(dist)

dist <- log(dist)

edist <- matrix(nrow = length(rownames(samdf2)), ncol = length(rownames(samdf2)))
edist

for(i in 1:nrow(edist)) for(j in 1:ncol(edist)) edist[i,j] <- samdf2$Ecotype[i] == samdf2$Ecotype[j]
boxplot(dist[lower.tri(dist)] ~ edist[lower.tri(edist)],
                    outline=FALSE,
                    ylab = "Bray-Curtis Pairwise Distance (ASVs)",
                    xlab = "Belong to Same Ecotype",
                    main = "Ecotype")

edist <- matrix(nrow = length(rownames(samdf2)), ncol = length(rownames(samdf2)))
edist

for(i in 1:nrow(edist)) for(j in 1:ncol(edist)) edist[i,j] <- samdf2$Ploid[i] == samdf2$Ploidy[j]
boxplot(dist[lower.tri(dist)] ~ edist[lower.tri(edist)],
                    outline=FALSE,
                    ylab = "Bray-Curtis Pairwise Distance (ASVs)",
                    xlab = "Belong to Same Ploidy",
                    main = "Ploidy")

edist <- matrix(nrow = length(rownames(samdf2)), ncol = length(rownames(samdf2)))
edist

for(i in 1:nrow(edist)) for(j in 1:ncol(edist)) edist[i,j] <- samdf2$Geo[i] == samdf2$Geo[j]
boxplot(dist[lower.tri(dist)] ~ edist[lower.tri(edist)],
                    outline=FALSE,
                    ylab = "Bray-Curtis Pairwise Distance (ASVs)",
                    xlab = "Belong to Same Genetic Cluster",
                    main = "Genetic Cluster")

# #### Generate Taxonomy Tables #### 
# # : https://github.com/joey711/phyloseq/issues/616
phyloGlom.c = tax_glom(ps.rare, "class")
glomTax.c = tax_table(phyloGlom.c)[,"class"]
glomOTU.c = otu_table(phyloGlom.c, taxa_are_rows=FALSE)
glomTable.c = merge(glomOTU.c,glomTax.c,by=0,all=TRUE)
rownames(glomTable.c) = glomTable.c[,"class"]
glomTable.c$Row.names = NULL
glomTable.c$class = NULL
class.tax <- phyloGlom.c@tax_table
rownames(class.tax) <- class.tax[,"class"]
ps.class <- phyloseq(otu_table(glomTable.c, taxa_are_rows=TRUE), tax_table(class.tax), sample_data(samdf))
ps.class
#
phyloGlom.p = tax_glom(ps.rare, "phylum")
glomTax.p = tax_table(phyloGlom.p)[,"phylum"]
glomOTU.p = otu_table(phyloGlom.p, taxa_are_rows=FALSE)
glomTable.p = merge(glomOTU.p,glomTax.p,by=0,all=TRUE)
rownames(glomTable.p) = glomTable.p[,"phylum"]
glomTable.p$Row.names = NULL
glomTable.p$phylum = NULL
phylum.tax <- phyloGlom.p@tax_table
rownames(phylum.tax) <- phylum.tax[,"phylum"]
ps.phylum <- phyloseq(otu_table(glomTable.p, taxa_are_rows=TRUE), tax_table(phylum.tax), sample_data(samdf))
ps.phylum
# 
phyloGlom.o = tax_glom(ps.rare, "order")
glomTax.o = tax_table(phyloGlom.o)[,"order"]
glomOTU.o = otu_table(phyloGlom.o, taxa_are_rows=FALSE)
glomTable.o = merge(glomOTU.o,glomTax.o,by=0,all=TRUE)
rownames(glomTable.o) = glomTable.o[,"order"]
glomTable.o$Row.names = NULL
glomTable.o$order = NULL
order.tax <- phyloGlom.o@tax_table
rownames(order.tax) <- order.tax[,"order"]
ps.order <- phyloseq(otu_table(glomTable.o, taxa_are_rows=TRUE), tax_table(order.tax), sample_data(samdf))
ps.order
# 
phyloGlom.f = tax_glom(ps.rare, "family")
glomTax.f = tax_table(phyloGlom.f)[,"family"]
glomOTU.f = otu_table(phyloGlom.f, taxa_are_rows=FALSE)
glomTable.f = merge(glomOTU.f,glomTax.f,by=0,all=TRUE)
rownames(glomTable.f) = glomTable.f[,"family"]
glomTable.f$Row.names = NULL
glomTable.f$family = NULL
family.tax <- phyloGlom.f@tax_table
rownames(family.tax) <- family.tax[,"family"]
ps.family <- phyloseq(otu_table(glomTable.f, taxa_are_rows=TRUE), tax_table(family.tax), sample_data(samdf))
fam_table <- t(ps.family@otu_table)
colnames(fam_table)
dim(fam_table)
write.csv(fam_table, "fam_table_365.csv")
#
phyloGlom.g = tax_glom(ps.rare, "genus")
glomTax.g = tax_table(phyloGlom.g)[,"genus"]
glomOTU.g = otu_table(phyloGlom.g, taxa_are_rows=FALSE)
glomTable.g = merge(glomOTU.g,glomTax.g,by=0,all=TRUE)
rownames(glomTable.g) = glomTable.g[,"genus"]
glomTable.g$Row.names = NULL
glomTable.g$genus = NULL
genus.tax <- phyloGlom.g@tax_table
rownames(genus.tax) <- genus.tax[,"genus"]
ps.genus <- phyloseq(otu_table(glomTable.g, taxa_are_rows=TRUE), tax_table(genus.tax), sample_data(samdf))


#### Plot Relative Abundances ####
glom <- tax_glom(ps.rare, taxrank = "phylum")
glom.merge <- merge_samples(glom, "Ecotype")
glom.merge.prop <- transform_sample_counts(glom.merge, function(x) x / sum(x))
glom.merge.prop.melt <- psmelt(glom.merge.prop) # create dataframe from phyloseq object
glom.merge.prop.melt.filt <- filter(glom.merge.prop.melt, Abundance > 0.02)
data_glom <- arrange(glom.merge.prop.melt.filt, phylum)
data_glom

data_glom$Sample <- factor(data_glom$Sample, levels = c("Lowland", "Both", "Upland"))
p <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill= phylum)) +
  geom_bar(stat="identity", position = "stack") +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition between Ecotypes") +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = phylum_colors) + theme(legend.position="bottom")
p

glom.merge <- merge_samples(glom, "Geo")
glom.merge.prop <- transform_sample_counts(glom.merge, function(x) x / sum(x))
glom.merge.prop.melt <- psmelt(glom.merge.prop) # create dataframe from phyloseq object
glom.merge.prop.melt.filt <- filter(glom.merge.prop.melt, Abundance > 0.02)
data_glom <- arrange(glom.merge.prop.melt.filt, phylum)
data_glom
data_glom$Sample <- factor(data_glom$Sample, levels = c("North", "South", "Northeast", "West", "East"))
p <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill= phylum)) +
  geom_bar(stat="identity", position = "stack") +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition between Genetic Clusters") +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = phylum_colors) + theme(legend.position="bottom")
p

glom.merge <- merge_samples(glom, "Ploidy")
glom.merge.prop <- transform_sample_counts(glom.merge, function(x) x / sum(x))
glom.merge.prop.melt <- psmelt(glom.merge.prop) # create dataframe from phyloseq object
glom.merge.prop.melt.filt <- filter(glom.merge.prop.melt, Abundance > 0.02)
data_glom <- arrange(glom.merge.prop.melt.filt, phylum)
data_glom
data_glom$Sample <- factor(data_glom$Sample, levels = c("Four", "Hybrid", "Eight"))
p <- ggplot(data=data_glom, aes(x=Sample, y=Abundance, fill= phylum)) +
  geom_bar(stat="identity", position = "stack") +
  ylab("Relative Abundance (Phyla > 2%) \n") +
  ggtitle("Phylum Composition between Ploidy Levels") +
  theme(axis.title.x = element_blank())+
  scale_fill_manual(values = phylum_colors) + theme(legend.position="bottom")
p

#### Core Microbiome Analysis ####
library(eulerr)
library(microbiome)

table(meta(ps.family)$Geo, useNA = "always")

pseq.rel <- microbiome::transform(ps.family, "compositional")

library(RColorBrewer)
prevalences <- seq(.05, 1, .05)
detections <- 10^seq(log10(1e-4), log10(.2), length = 10)

p1 <- plot_core(pseq.rel, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(5, "RdBu")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")

p1 <- p1 + theme_bw() + ylab("ASVs")
p1 + theme(axis.text.x = element_text(angle = 45))


core90 <- core(pseq.rel, detection = 0, prevalence = 0.9)
core10 <- core(pseq.rel, detection = 0, prevalence = 0.1)

taxa(core90)
taxa(core10)

fams <- rownames(pseq.rel@tax_table)

prevs <- as.data.frame(fams)
prevs$core90 <- prevs$fam %in% taxa(core90)
prevs$core10 <- prevs$fam %in% taxa(core10)

write.csv(prevs, 'prevs.csv')

t(ps.family@otu_table)[1:5,1:5]

write.csv(t(ps.family@otu_table), "family_otus.csv")

#### Venn Diagrams ####

#Genetic Cluster
library(ggpubr)
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

mycols <- c(West = "#ABDDDE", South = "#096E9C", East = "#2E3436", North = "#C93312", Northeast = "#D79C4E") 

venn_geo <- plot(venn(list_core),
     fills = mycols,
     main = "Core Microbiome - Genetic Cluster", cex = 4) 


#Ecotype

eco <- unique(as.character(meta(pseq.rel)$Ecotype))
print(eco)

list_core <- c() # an empty object to store information

for (n in eco){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Ecotype == n) # Choose sample from Ecotype by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, # 0.001 in atleast 90% samples 
                         prevalence = 0.9)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)

mycols <- c(phylum_colors[31],phylum_colors[30],phylum_colors[7]) 
venn_eco <- plot(venn(list_core),
     fills = mycols,
     main = "Core Microbiome - Ecotype")

#Ploidy

ploidy <- unique(as.character(meta(pseq.rel)$Ploidy))
print(ploidy)

list_core <- c() # an empty object to store information

for (n in ploidy){ # for each variable n in DiseaseState
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq.rel, Ploidy == n) # Choose sample from Ecotype by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0, # 0.001 in atleast 90% samples 
                         prevalence = 0.9)
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

print(list_core)
mycols <- c(phylum_colors[16],phylum_colors[18],phylum_colors[21]) 

mycols <- c(phylum_colors[21],phylum_colors[16],phylum_colors[18]) 
venn_ploidy<- plot(venn(list_core),
     fills = mycols,
     main = "Core Microbiome - Ploidy")


library("cowplot")
ggarrange(venn_eco, venn_ploidy, venn_geo, labels = c("A", "B", "C"), ncol =3)



