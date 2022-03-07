#### Jesse's Server with GWASanth /.RData
setwd("~/switchgrass/gwas_core_taxa")
load("~/switchgrass/gwas_core_taxa/.RData")

library(statgenGWAS)
library(vroom)
library(textshape)
library(dplyr)

tet <- vroom("tetra_Gnum_MAF.csv", col_names = T)
tet[1:5,1:5]
tet <- textshape::column_to_rownames(tet, loc = 1)
fam <- read.csv("family_otus.csv", header = T, row.names = 1)

setwd("~/switchgrass/gwas_core_taxa/fours_only")

# tet <- tet2
tet2 <- tet
# write.csv(tet2, "ALL_Gnum_maf_361.csv")
map <- map2
head(map)
colnames(map) <- c("chr", "pos")
map$SNP.names <- row.names(map)
map <- map[, c(3,1,2)]
map$SNP.names <- sub("^", "SNP", map$SNP.names )
rownames(map) <- map$SNP.names
head(map)
dim(map)

fournames <- c("Actinomycetaceae", "Alcanivoracaceae1", "Amb.16S.1034", "Aneurinibacillaceae", "Barnesiellaceae", "Bogoriellaceae", "Brevibacillaceae", "Caedibacteraceae", "Caldilineaceae", "Coxiellaceae", "Dermabacteraceae", "Desulfallas.Sporotomaculum", "Desulfitobacteriaceae", "Desulfobulbaceae", "Desulfovibrionaceae", "DEV007", "Hydrogenophilaceae", "Idiomarinaceae", "Kaistiaceae", "LiUU.11.161", "Magnetospirillaceae", "Methyloligellaceae", "Neisseriaceae", "Oligosphaeraceae", "Oxobacteraceae", "Paludibacteraceae", "Porphyromonadaceae", "Pseudanabaenaceae", "Psychromonadaceae", "Puniceicoccaceae", "Saccharimonadaceae", "Sanguibacteraceae", "Selenomonadaceae", "Shewanellaceae", "Sumerlaeaceae", "Tannerellaceae", "Thermaceae", "Thermomicrobiaceae", "Trueperaceae", "UCG.010", "Wohlfahrtiimonadaceae")
dim(fam)
names(fam)

fam <- fam[,fournames]
fam <- na.omit(fam)
dim(fam)
fam <- cbind(genotype = rownames(fam), fam)

tnames <- colnames(tet)
head(tnames)
tnames <- sub("X", "", tnames)
tnames <- as.numeric(tnames)
SNPnames <- sub("^", "SNP", tnames)
colnames(tet) <- SNPnames
tet[1:5,1:5]

fam[1:5,1:5]
tet[1:5,1:5]
map[1:5,1:3]

dim(map)
sapply(map, class)
map[2] <- lapply(map[2], as.character)
map[3] <- lapply(map[3], as.numeric)

sapply(map, class)

sapply(tet[1:5,1:5], class)
sapply(fam, class)
dim(fam)
fam[2:42] <- lapply(fam[2:42], as.numeric)
fam[3] <- lapply(fam[3], as.numeric)

tet <- tet[ order(row.names(tet)), ]

tet[1:5,1:5]
tet[tet == 1] <- 3
tet[tet == 0] <- 1
tet[tet == -1] <- 0
tet[tet == 3] <- 2
tet[1:5,1:5]


tet <- t(tet)
tet <- subset(tet, rownames(tet) %in% rownames(map))
map <- subset(map, rownames(map) %in% rownames(tet))
tet <- t(tet)
tet[1:5,1:5]
dim(tet)
dim(map)

tet <- subset(tet, rownames(tet) %in% rownames(fam))
fam <- subset(fam, rownames(fam) %in% rownames(tet))


dim(map)
dim(tet)
dim(fam)

fam2 <- fam
fam2 <-  fam2[, colSums(fam2 != 0) > 1]
dim(fam)
dim(fam2)
fam2
gDataFam <- createGData(geno = tet, map = map, pheno = fam2)
set.seed(1)
gDataFamDedup <- codeMarkers(gDataFam, impute = FALSE, verbose = TRUE)

rm(GWASfam)
GWASfam <- runSingleTraitGwas(gData = gDataFam,
                              thrType = "fdr",
                              GLSMethod = "multi",
                              sizeInclRegion = 25000,
                              pThr = 0.01,
                              kinshipMethod = "vanRaden",
                              MAF = 0.05)

summary(GWASfam)

print(GWASfam$signSnp, row.names = T)

library(ggplot2)
plot(GWASfam, plotType = "manhattan", trait = "Xanthobacteraceae") + ggtitle("Xanthobacteraceae")
plot(GWASfam, plotType = "manhattan", trait = "Sphingomonadaceae") + ggtitle("Sphingomonadaceae")
plot(GWASfam, plotType = "manhattan", trait = "Comamonadaceae") + ggtitle("Comamonadaceae")
plot(GWASfam, plotType = "manhattan", trait = "Rhodanobacteraceae") + ggtitle("Rhodanobacteraceae")
plot(GWASfam, plotType = "manhattan", trait = "Mycobacteriaceae") + ggtitle("Mycobacteriaceae")
plot(GWASfam, plotType = "manhattan", trait = "Nitrosomonadaceae") + ggtitle("Nitrosomonadaceae")
plot(GWASfam, plotType = "manhattan", trait = "Chthoniobacteraceae") + ggtitle("Chthoniobacteraceae")
plot(GWASfam, plotType = "manhattan", trait = "Gaiellaceae") + ggtitle("Gaiellaceae")
plot(GWASfam, plotType = "manhattan", trait = "Nitrospiraceae") + ggtitle("Nitrospiraceae")
plot(GWASfam, plotType = "manhattan", trait = "Pedosphaeraceae") + ggtitle("Pedosphaeraceae")
plot(GWASfam, plotType = "manhattan", trait = "Comamonadaceae") + ggtitle("Comamonadaceae")
plot(GWASfam, plotType = "manhattan", trait = "Chitinophagaceae") + ggtitle("Chitinophagaceae")
plot(GWASfam, plotType = "manhattan", trait = "Pyrinomonadaceae") + ggtitle("Pyrinomonadaceae") 


plot(GWASfam, plotType = "qq", trait = "Xanthobacteraceae") + ggtitle("Xanthobacteraceae") 
plot(GWASfam, plotType = "qq", trait = "Chitinophagaceae") + ggtitle("Chitinophagaceae")
plot(GWASfam, plotType = "qq", trait = "Sphingomonadaceae") + ggtitle("Sphingomonadaceae")
plot(GWASfam, plotType = "qq", trait = "Pyrinomonadaceae") + ggtitle("Pyrinomonadaceae") 

write.csv(GWASfam$signSnp, "core_snps_tet_all_fam.csv")


#### Average Anth ####

phenoList2 <- split(x = pheno[c("genotype", "Height_cm",
                                "Circumference_cm", "Volume_cm3",
                                "Anthracnose", "Vigor")], 
                    f = pheno[c("Year")])

head(pheno)
head(phenoList2)

mean_2019 <- aggregate(phenoList2$`2019`[, 5], list(phenoList2$`2019`$genotype), mean)

mean_2019

colnames(mean_2019) <- c("genotype", "Anth_avg")
head(mean_2019)
mean_2019$ALL <- "average"
head(mean_2019)

mean_2020 <- aggregate(phenoList2$`2020`[, 5], list(phenoList2$`2020`$genotype), mean)

mean_2020

colnames(mean_2020) <- c("genotype", "Anth_avg")
head(mean_2020)
mean_2020$ALL <- "average"
head(mean_2020)


pheno.mean.2019 <- split(x = mean_2019[c("genotype", "Anth_avg")], 
                         f = mean_2019[["ALL"]])

pheno.mean.2020 <- split(x = mean_2020[c("genotype", "Anth_avg")], 
                         f = mean_2019[["ALL"]])

gData2019 <- createGData(tet = tet, map = map, pheno = pheno.mean.2019)
gData2020 <- createGData(tet = tet, map = map, pheno = pheno.mean.2020)


gData2019Dedup <- codeMarkers(gData2019, impute = FALSE, verbose = TRUE)
gData2020Dedup <- codeMarkers(gData2020, impute = FALSE, verbose = TRUE)


str(gData2Dedup)

GWAS2019 <- runSingleTraitGwas(gData = gData2019Dedup,
                               trials = "average",
                               traits = "Anth_avg",
                               thrType = "fdr",
                               GLSMethod = "multi",
                               pThr = 0.1,
                               kinshipMethod = "vanRaden")
summary(GWAS2019)

GWAS2020 <- runSingleTraitGwas(gData = gData2020Dedup,
                               trials = "average",
                               traits = "Anth_avg",
                               thrType = "fdr",
                               GLSMethod = "multi",
                               pThr = 0.1,
                               kinshipMethod = "vanRaden")
summary(GWAS2020)

plot(GWAS2020, plotType = "manhattan", trial = "average", trait = "Anth_avg") + ggtitle("Avg. Anthracnose 2020")
plot(GWAS2019, plotType = "manhattan", trial = "average", trait = "Anth_avg") + ggtitle("Avg. Anthracnose 2019")
plot(GWAS2, plotType = "manhattan", trial = "average", trait = "Anth_avg") + ggtitle("Avg. Anthracnose 2019-2020")


#### GO Analysis ####

# GO <- read.csv("GO_terms_tet_LD.csv")
GO

library(ggplot2)
#"GO.term"              "Ontology"             "Description"          "Number.in.input.list"
#"Number.in.BG.Ref"     "p.value"              "FDR"

GO$perc <- GO$Number.in.input.list / sum(GO$Number.in.input.list)
GO$perc <- GO$perc * 100

GO$perc2 <- GO$Number.in.BG.Ref / sum(GO$Number.in.BG.Ref)
GO$perc2 <- GO$perc2 * 100

ont <- GO$Ontology

ont <- gsub("C", "Cellular Component", ont)
ont <- gsub("P", "Biological Process", ont)
ont <- gsub("F", "Molecular Function", ont)
ont
GO$Ontology <- ont
GO$Ontology <- factor(GO$Ontology, levels = c("Molecular Function", "Biological Process", "Cellular Component"))

GOplot <- ggplot(GO, aes(x= reorder(Description, -perc), y=perc, fill=Ontology)) +
  geom_bar(stat="identity") + theme_light() + scale_y_continuous(labels=function(x) paste0(x,"%"))

GOplot + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7)) +
  ggtitle("GO Terms") +
  xlab("GO Term Description") + ylab("Percent of Genes") +
  scale_fill_manual(values=c("#09397e","#655e03","#D4A05B"))



#### Sig SNPS ####

sigs <- read.csv("sig_snps_tet_processed.csv")
head(sigs)

sigs_fam <- sigs['fam.trait']
head(sigs_fam)

sig_table <- as.data.frame(table(sigs_fam))
sig_table <- textshape::column_to_rownames(sig_table, loc = 1)
sig_table <- as.matrix(sig_table)


h2 <- read.csv("h2_fam.csv")
h2 <- textshape::column_to_rownames(h2, loc = 1)
h2 <- as.matrix(h2)


new_names <- intersect(rownames(sig_table), rownames(h2))
h2 <- subset(h2, rownames(h2) %in% new_names)
sig_table <- subset(sig_table, rownames(sig_table) %in% new_names)
h2 <- h2[ order(row.names(h2)), ]
sig_table <- sig_table[ order(row.names(sig_table)), ]



sig_h2 <- cbind(sig_table, h2)
tail(sig_h2)

sig_h2 <- as.data.frame(sig_h2)
head(sig_h2)

prev90 <- c("Xanthobacteraceae", "Sphingomonadaceae", "Comamonadaceae", "Rhodanobacteraceae", "Nitrospiraceae", 
            "Haliangiaceae", "Mycobacteriaceae", "Nitrosomonadaceae", "Chthoniobacteraceae", "Pedosphaeraceae",
            "Chitinophagaceae", "Gaiellaceae", "Pyrinomonadaceae")

prev30 <- c("Xiphinematobacteraceae", "Xanthomonadaceae", "Xanthobacteraceae", "WD2101,soil,group", "Vicinamibacteraceae", "Vermiphilaceae", "Unknown,Family,4", "TRA3-20", "Thermoanaerobaculaceae", "Streptomycetaceae", "Steroidobacteraceae", "Sphingomonadaceae", "Sphingobacteriaceae", "Solirubrobacteraceae", "Solibacteraceae", "SM2D12", "SC-I-84", "Rhodanobacteraceae", "Rhizobiales,Incertae,Sedis", "Rhizobiaceae", "Reyranellaceae", "Pyrinomonadaceae", "Pseudonocardiaceae", "Pseudomonadaceae", "Polyangiaceae", "Pirellulaceae", "Phycisphaeraceae", "Phaselicystidaceae", "Pedosphaeraceae", "Oxalobacteraceae", "Opitutaceae", "Nocardioidaceae", "Nitrospiraceae", "Nitrosomonadaceae", "Nakamurellaceae", "Myxococcaceae", "Mycobacteriaceae", "Microscillaceae", "Micromonosporaceae", "Micrococcaceae", "Microbacteriaceae", "Methylophilaceae", "Labraceae", "Ktedonobacteraceae", "Koribacteraceae", "Kineosporiaceae", "KF-JG30-B3", "JG30-KF-AS9", "Isosphaeraceae", "Intrasporangiaceae", "Iamiaceae", "Hyphomicrobiaceae", "Haliangiaceae", "Geodermatophilaceae", "Gemmatimonadaceae", "Gemmataceae", "Geminicoccaceae", "Gaiellaceae", "Frankiaceae", "Flavobacteriaceae", "env.OPS,17", "Dongiaceae", "Diplorickettsiaceae", "Devosiaceae", "cvE6", "Comamonadaceae", "Clostridiaceae", "Chthonomonadaceae", "Chthoniobacteraceae", "Chitinophagaceae", "Caulobacteraceae", "Burkholderiaceae", "Bryobacteraceae", "Blastocatellaceae", "BIrii41", "Beijerinckiaceae", "Bdellovibrionaceae", "Bacillaceae", "Azospirillaceae", "Anaeromyxobacteraceae", "Alicyclobacillaceae", "AKYH767", "Acidothermaceae", "Acidobacteriaceae,(Subgroup,1)", "Acetobacteraceae", "A21b", "X67-14")

sig_h2 <- cbind(genotype = rownames(sig_h2), sig_h2)

sig_h2$prev90 <- as.numeric(sig_h2$genotype %in% prev90)
sig_h2$prev30 <- as.numeric(sig_h2$genotype %in% prev30)

head(sig_h2)

final <- sig_h2 %>% select(sig_table, h2, se, prev30, prev90)
final <- final %>% rename(SNP_count = sig_table)
head(final)
write.csv(final, "SNPcount_h2_prev.csv")
