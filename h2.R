setwd("~/switchgrass/heritability")
load("~/switchgrass/heritability/.RData")
library(sommer)
library(vroom)
library(textshape)
library(dplyr)
library(magicfor)  
library(splitstackshape)
options(expressions = 5e5)


Gnum <- vroom("ALL_Gnum_MAF_362_processed.csv", col_names = T)
y <- read.csv("fam_table_365.csv", row.names = 1, header = T)

Gnum[1:25,1:5]
Gnum <- textshape::column_to_rownames(Gnum, loc = 1)
Gnum <- Gnum[ order(row.names(Gnum)), ]
Gnum[1:5,1:5]

y[1:5,1:5]
dim(y)
dim(Gnum)

rownames(y)
rownames(Gnum)
 
Gnum = Gnum[row.names(Gnum) %in% row.names(y), ]
y = y[row.names(y) %in% row.names(Gnum), ]
dim(y)
dim(Gnum)

y <- y[ order(row.names(y)),]

Gnum[1:5,1:5]
y[1:5,1:5]


#Kinship Matrices
rm(A,D)
Gnum <- as.matrix(Gnum)
A <- A.mat(Gnum) #additive
D <- D.mat(Gnum) #dominant (intra-locus interactions)
# E <- E.mat(Gnum) #epistatic (inter-locus interactions) -  Ignore this one


famnames <- colnames(y)

y$id <- row.names(y)
y$idd <- y$id
y$ide <- y$id


fam <- y
magic_for(print, silent = TRUE) # call magic_for()
for(i in famnames){
  print(i)
  fam$tmp <- fam[,i]
  ans.ADE.x <-mmer(tmp ~ 1 , random = ~ vs(id,Gu=A) + vs(idd,Gu=D), rcov=~units, tolparinv = 10, data=fam, verbose = F)
  print(vpredict(ans.ADE.x, h2~(V1)/( V1+V3)))
}

out <- magic_result_as_dataframe()
out <- cSplit(out, "vpredict(ans.ADE.x,h2~(V1)/(V1+V3))", sep=",")
out <- as.data.frame(out)
names(out) <- c("i", "taxa", "h2", "se")
out$h2 <- gsub(".*= ", "", out$h2)
out$se <- gsub(".*= ", "", out$se)
out$se <- gsub(")", "", out$se)
out$i <- NULL
out$h2 <- as.numeric(out$h2)
out$se <- as.numeric(out$se)
out$tf <- with(out, h2 - se > 0)
out$level <- paste("ASV")
out$site <- paste("ALL")
out <- out[order(-out$tf),] 
head(out)

write.csv(out,"h2_fam.csv", row.names = FALSE)

#### Plot h2 ####
# output from above filtered in Excel where (h2 - se) > 0 = "core_h2_scores.csv" 

library(ggplot2)
library(ggrepel)
theme_set(theme_light())

set.seed(42)

h2 <- read.csv('core_h2_scores.csv', header = T)

head(h2)

plot <- ggplot(h2, aes(x=reorder(family, -h2), y=h2)) + 
  geom_errorbar(aes(ymin = h2 - se, ymax = h2 + se), width=0.2) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(10, 25, 10, 10)) + 
  geom_point() +
  ylab(expression(paste(italic("h")^"2"))) + 
  xlab("Family") +
  ggtitle("Heritable Core Families") + 
  scale_y_continuous(labels = scales::percent) 


plot + theme(axis.text.x = element_text(color = h2$color))



             