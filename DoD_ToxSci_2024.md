DoD_single_vs_repeated
================
Arun Ghosh
2024-11-20

``` r
# -loading packages and data

library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(dplyr)
library(stringr)
library(ggplot2)
library(EnhancedVolcano)
library(eulerr)
library(purrr)
library(ggVennDiagram)
library(tidyverse)
library(pheatmap)
library(RColorBrewer) 
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(clusterProfiler)
library(readxl)
library(psych)
library(corrplot)
library(enrichplot)
library(openxlsx)
library(DESeq2)
library(data.table)
library(quanteda)
library(eulerr)
library(DOSE)
```

``` r
counts <- data.frame(read_excel("raw_counts.xlsx", sheet = "raw_counts"))

colnames(counts)[1] <- "Geneid"
```

``` r
##  Donor effect correction and analysis #########

d0 <- DGEList(counts)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 10 # genes expressed in at least 10 samples to be included
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(counts) # Sample names
snames <- snames[2:length(snames)]

group <- substr(snames, 10, (nchar(snames))) # for exposure groups
group <- as.factor(group)
group

donor <- substr(snames, 6, 8) # for donor based batch correction
batch <- as.factor(donor)
batch

mm <- model.matrix(~0+group+batch)

y <- voom(d, mm, plot = T)
```

![](README_figs/README--Limma%20Voom-1.png)<!-- -->

``` r
colnames(mm)
fit <- lmFit(y, mm) # Fitting linear models in limma

head(coef(fit))

x <- colnames(coef(fit))
length(x)

x # to see the groups

#------------------------------------------------------------#
```

``` r
#------------------------------------------------------------#
# for comparison within single exposure groups

x1 <- char_select(x, "*Mock_1", valuetype = "glob") # selecting levels
# to use in the "for" loop

a1 <- list() # list of both coding and non-coding transcripts
b1 <- list() # for storing analyzed data fro significantly altered genes
c1 <- list() # for getting the list of all genes
d1 <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x1)){if(x1[i] != "groupPBS_Mock_1"){
  difference <- paste(x1[i],"-","groupPBS_Mock_1", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x1[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  c1[[i]] <- top.table
  names(c1)[i] <- exposure
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d1[[i]] = top.table$Geneid
  names(d1)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a1[[i]] = top.table$symbol
  names(a1)[i] <- exposure
  b1[[i]] <- top.table
  names(b1)[i] <- exposure}}


# removing the empty control group from the lists
a1[3] <- NULL
b1[3] <- NULL
c1[3] <- NULL
d1[3] <- NULL
```

``` r
#------------------------------------------------------------#
# for comparison within repeated exposure groups

x2 <- char_select(x, "*Mock_3", valuetype = "glob") # selecting levels
# to use in the "for" loop

a2 <- list() # list of both coding and non-coding transcripts
b2 <- list() # for storing analyzed data for significantly altered genes
c2 <- list() # for getting the list of all genes
d2 <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x2)){if(x2[i] != "groupPBS_Mock_3"){
  difference <- paste(x2[i],"-","groupPBS_Mock_3", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x2[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  c2[[i]] <- top.table
  names(c2)[i] <- exposure
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d2[[i]] = top.table$Geneid
  names(d2)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a2[[i]] = top.table$symbol
  names(a2)[i] <- exposure
  b2[[i]] <- top.table
  names(b2)[i] <- exposure}}


# removing the empty control group from the lists
a2[3] <- NULL
b2[3] <- NULL
c2[3] <- NULL
d2[3] <- NULL
```

``` r
EP_1X <- euler(a1, shape = "ellipse") # Euler plot 

plot(EP_1X,
     quantities = TRUE, cex = 5,
     main = "single",
     labels = list(font = 3, cex = 1),
     fills = c("lightgreen", "white","lightblue", "white")
)
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-1.png)<!-- -->

``` r
EP_3X <- euler(a2, shape = "circle") # Euler plot 
plot(EP_3X,
     quantities = TRUE, cex = 5,
     main = "repeated",
     labels = list(font = 3, cex = 1),
     fills = c("lightgreen", "white","lightblue", "white")
)
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-2.png)<!-- -->

``` r
#################################################################
# comparison between single and repeated exposures

p <- list()
p[[1]]<- a1$CF
names(p)[1] <- "single"

p[[2]]<- a2$CF
names(p)[2] <- "repeated"

EP_CF <- euler(p, shape = "circle") # Euler plot 
plot(EP_CF,
     quantities = TRUE, cex = 2,
     main = "Cardboard Flaming -single vs repeated",
     lty = 1:1,
     labels = list(font = 3, cex = 2),
     fills = c( "white", "aquamarine" ))# "burlywood1"
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-3.png)<!-- -->

``` r
############################################################
p <- list()
p[[1]]<- a1$CS
names(p)[1] <- "single"

p[[2]]<- a2$CS
names(p)[2] <- "repeated"

p <- ggVennDiagram(p, label = "count") + ggtitle("Cardboard Smoldering")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-4.png)<!-- -->

``` r
############################################################
p <- list()
p[[1]]<- a1$PF
names(p)[1] <- "single"

p[[2]]<- a2$PF
names(p)[2] <- "repeated"

EP_PF <- euler(p, shape = "circle") # Euler plot 
plot(EP_PF,
     quantities = TRUE, cex = 2,
     main = "Plastic Flaming -single vs repeated",
     lty = 1:1,
     labels = list(font = 3, cex = 2),
     fills = c( "white", "aquamarine" ))
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-5.png)<!-- -->

``` r
############################################################
p <- list()
p[[1]]<- a1$PS
names(p)[1] <- "single"

p[[2]]<- a2$PS
names(p)[2] <- "repeated"


p <- ggVennDiagram(p, label = "count") + ggtitle("Plastic Smoldering")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams-6.png)<!-- -->

``` r
##########################################################################

Gdata <- c2 # repeated

#---------------------------------------------#
# Cardboard exposure group -Flaming

Gdata$CF <- Gdata$CF %>% mutate(ProbeID = Gdata$CF$Geneid)   

GdataCBf <- Gdata$CF #Cardboard
GdataCBf <- select(GdataCBf, ProbeID, logFC)
rownames(GdataCBf) <- NULL

#making ranked gene list
genelist_GdataCBf = GdataCBf[,2] #numeric vector
names(genelist_GdataCBf) = as.character(GdataCBf[,1]) #named vector
genelist_GdataCBf = sort(genelist_GdataCBf, decreasing = TRUE) #must sort in descending order


#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#

#GO over-representation analysis
geneCBf <- names(genelist_GdataCBf[genelist_GdataCBf[] > 0.5]) # logFC > 0.5
enGOCBf <- enrichGO(gene         = geneCBf,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)

#head(enGOCBf)
enGOCBf_c3 <- gsfilter(enGOCBf, by = 'Count', min = 3) # selecting categories with at least gene count 3
#barplot
barplot(enGOCBf_c3, showCategory=10, font = 15, title = "Cardboard Flaming -repeated")+ 
  scale_x_continuous(breaks = seq(0, 10, by = 5), limits=c(0,10))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold")) +
  scale_y_discrete(label=function(y) stringr::str_trunc(y, 40))
```

![](README_figs/README--GO%20enrichment-1.png)<!-- -->

``` r
##################################################################
#---------------------------------------------#
# Plastic exposure group -Flaming

Gdata$PF <- Gdata$PF %>% mutate(ProbeID = Gdata$PF$Geneid)   

GdataPLf <- Gdata$PF #Plastic
GdataPLf <- select(GdataPLf, ProbeID, logFC)
rownames(GdataPLf) <- NULL

#making ranked gene list
genelist_GdataPLf = GdataPLf[,2] #numeric vector
names(genelist_GdataPLf) = as.character(GdataPLf[,1]) #named vector
genelist_GdataPLf = sort(genelist_GdataPLf, decreasing = TRUE) #must sort in descending order


#--------------------------------------------------------------------#
#--------------------------------------------------------------------#

#GO over-representation analysis
genePLf <- names(genelist_GdataPLf[genelist_GdataPLf[] > 0.5]) # logFC > 0.5
enGOPLf <- enrichGO(gene         = genePLf,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)

#head(enGOPLf)
enGOPLf_c3 <- gsfilter(enGOPLf, by = 'Count', min = 3) # selecting categories with at least gene count 3
#barplot
barplot(enGOPLf_c3, showCategory=10, font = 15, title = "Plastic Flaming -repeated")+ 
  scale_x_continuous(breaks = seq(0, 10, by = 5), limits=c(0,10))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))+
  scale_y_discrete(label=function(y) stringr::str_trunc(y, 40))
```

![](README_figs/README--GO%20enrichment-2.png)<!-- -->

``` r
##########################################################################

Gdata <- c1 # single

# #---------------------------------------------#
# # Cardboard exposure group -Flaming [NO RESULT]
######################################-----------#

# Gdata$CF <- Gdata$CF %>% mutate(ProbeID = Gdata$CF$Geneid)   
# 
# GdataCBf <- Gdata$CF #Cardboard
# GdataCBf <- select(GdataCBf, ProbeID, logFC)
# rownames(GdataCBf) <- NULL
# 
# #making ranked gene list
# genelist_GdataCBf = GdataCBf[,2] #numeric vector
# names(genelist_GdataCBf) = as.character(GdataCBf[,1]) #named vector
# genelist_GdataCBf = sort(genelist_GdataCBf, decreasing = TRUE) #must sort in descending order
# 
# #--------------------------------------------------------------------#
# 
# #GO over-representation analysis
# geneCBf <- names(genelist_GdataCBf[genelist_GdataCBf[] > 0.5]) # logFC > 0.5
# enGOCBf <- enrichGO(gene         = geneCBf,
#                     OrgDb         = org.Hs.eg.db,
#                     keyType       = 'ENSEMBL',
#                     ont           = "ALL",
#                     pAdjustMethod = "BH",
#                     pvalueCutoff  = 0.05,
#                     qvalueCutoff  = 0.1)
# 
# 
# enGOCBf_c3 <- gsfilter(enGOCBf, by = 'Count', min = 3) # selecting categories with at least gene count 3
# 
# #barplot
# barplot(enGOCBf_c3, showCategory=10, font = 15, title = "Cardboard Flaming -1X")+ 
#   scale_x_continuous(breaks = seq(0, 3, by = 1), limits=c(0,3))+
#   theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
#         title = element_text(size = 15, face="bold"))+
#   scale_y_discrete(label=function(y) stringr::str_trunc(y, 35))


##################################################################
#---------------------------------------------#
# Plastic exposure group -Flaming

Gdata$PF <- Gdata$PF %>% mutate(ProbeID = Gdata$PF$Geneid)   

GdataPLf <- Gdata$PF #Plastic
GdataPLf <- select(GdataPLf, ProbeID, logFC)
rownames(GdataPLf) <- NULL

#making ranked gene list
genelist_GdataPLf = GdataPLf[,2] #numeric vector
names(genelist_GdataPLf) = as.character(GdataPLf[,1]) #named vector
genelist_GdataPLf = sort(genelist_GdataPLf, decreasing = TRUE) #must sort in descending order


#--------------------------------------------------------------------#

#GO over-representation analysis
genePLf <- names(genelist_GdataPLf[genelist_GdataPLf[] > 0.5]) # logFC > 0.5
enGOPLf <- enrichGO(gene         = genePLf,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)


#barplot


enGOPLf_c3 <- gsfilter(enGOPLf, by = 'Count', min = 3) # selecting categories with at least gene count 3

barplot(enGOPLf_c3, showCategory=10, font = 15, title = "Plastic Flaming -single")+ 
  scale_x_continuous(breaks = seq(0, 10, by = 5), limits=c(0,10))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))+
  scale_y_discrete(label=function(y) stringr::str_trunc(y, 40))
```

![](README_figs/README--GO%20enrichment-3.png)<!-- -->

``` r
##  Donor effect correction and analysis of FEMALE DONORS only #######
###########################################===============############

countsF <- counts


countsF <- countsF %>% select(contains("Geneid") | 
                                contains("15R") |
                                contains("40R") |
                                contains("57N") )

d0 <- DGEList(countsF)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 10 # genes expressed in at least 10 samples to be included
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(countsF) # Sample names
snames <- snames[2:length(snames)]

donor <- substr(snames, 6, 8) # for donor based batch correction
batch <- as.factor(donor)
batch

group <- substr(snames, 10, (nchar(snames))) # for exposure groups
#group <- interaction(sex, group)
group <- as.factor(group)
group 

#---------------------------------------------------------#

mm <- model.matrix(~0+group+batch)

y <- voom(d, mm, plot = T)
```

![](README_figs/README--FEMALE%20DONORS-1.png)<!-- -->

``` r
colnames(mm)
fit <- lmFit(y, mm) # Fitting linear models in limma

head(coef(fit))

x <- colnames(coef(fit))
length(x)

x # to see the groups

#------------------------------------------------------------#
#------------------------------------------------------------#
# for comparison within single

x1 <- char_select(x, "*Mock_1", valuetype = "glob") # selecting levels
# to use in the "for" loop

a3Fa <- list() # list of both coding and non-coding transcripts
b3Fa <- list() # for storing analyzed data fro significantly altered genes
d3Fa <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x1)){if(x1[i] != "groupPBS_Mock_1"){
  difference <- paste(x1[i],"-","groupPBS_Mock_1", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x1[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d3Fa[[i]] = top.table$Geneid
  names(d3Fa)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a3Fa[[i]] = top.table$symbol
  names(a3Fa)[i] <- exposure
  b3Fa[[i]] <- top.table
  names(b3Fa)[i] <- exposure}}


# removing the empty control group from the lists
a3Fa[3] <- NULL
b3Fa[3] <- NULL
d3Fa[3] <- NULL

#------------------------------------------------------------#
# for comparison within repeated

x2 <- char_select(x, "*Mock_3", valuetype = "glob") # selecting levels
# to use in the "for" loop

a3Fb <- list() # list of both coding and non-coding transcripts
b3Fb <- list() # for storing analyzed data for significantly altered genes
c3Fb <- list() # for getting the list of all genes
d3Fb <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x2)){if(x2[i] != "groupPBS_Mock_3"){
  difference <- paste(x2[i],"-","groupPBS_Mock_3", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x2[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  c3Fb[[i]] <- top.table
  names(c3Fb)[i] <- exposure
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d3Fb[[i]] = top.table$Geneid
  names(d3Fb)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a3Fb[[i]] = top.table$symbol
  names(a3Fb)[i] <- exposure
  b3Fb[[i]] <- top.table
  names(b3Fb)[i] <- exposure}}


# removing the empty control group from the lists
a3Fb[3] <- NULL
b3Fb[3] <- NULL
c3Fb[3] <- NULL
d3Fb[3] <- NULL
```

``` r
##  Donor effect correction and analysis of MALE DONORS only #######
###########################################=============##############

countsM <- counts

countsM <- countsM %>% select(contains("Geneid") | 
                                contains("68J") |
                                contains("21O") |
                                contains("67L") )

d0 <- DGEList(countsM)# Create DGEList object
d0 <- calcNormFactors(d0)
cutoff <- 10 # genes expressed in at least 10 samples to be included
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left

snames <- colnames(countsM) # Sample names
snames <- snames[2:length(snames)]

donor <- substr(snames, 6, 8) # for donor based batch correction
batch <- as.factor(donor)
batch

group <- substr(snames, 10, (nchar(snames))) # for exposure groups
#group <- interaction(sex, group)
group <- as.factor(group)
group 

#---------------------------------------------------------#

mm <- model.matrix(~0+group+batch)

y <- voom(d, mm, plot = T)
```

![](README_figs/README--MALE%20DONORS-1.png)<!-- -->

``` r
colnames(mm)
fit <- lmFit(y, mm) # Fitting linear models in limma

head(coef(fit))

x <- colnames(coef(fit))
length(x)

x # to see the groups

#------------------------------------------------------------#
#------------------------------------------------------------#
# for comparison within single

x1 <- char_select(x, "*Mock_1", valuetype = "glob") # selecting levels
# to use in the "for" loop

a3Ma <- list() # list of both coding and non-coding transcripts
b3Ma <- list() # for storing analyzed data fro significantly altered genes
d3Ma <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x1)){if(x1[i] != "groupPBS_Mock_1"){
  difference <- paste(x1[i],"-","groupPBS_Mock_1", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x1[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d3Ma[[i]] = top.table$Geneid
  names(d3Ma)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a3Ma[[i]] = top.table$symbol
  names(a3Ma)[i] <- exposure
  b3Ma[[i]] <- top.table
  names(b3Ma)[i] <- exposure}}


# removing the empty control group from the lists
a3Ma[3] <- NULL
b3Ma[3] <- NULL
d3Ma[3] <- NULL

#------------------------------------------------------------#
# for comparison within repeated

x2 <- char_select(x, "*Mock_3", valuetype = "glob") # selecting levels
# to use in the "for" loop

a3Mb <- list() # list of both coding and non-coding transcripts
b3Mb <- list() # for storing analyzed data for significantly altered genes
d3Mb <- list() # for storing ENSEMBL transcript names

for(i in 1:length(x2)){if(x2[i] != "groupPBS_Mock_3"){
  difference <- paste(x2[i],"-","groupPBS_Mock_3", sep="")
  contr <- makeContrasts(difference, levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  tmp <- eBayes(tmp)
  top.table <- topTable(tmp, sort.by = "P", n = Inf)
  exposure <- substr(x2[i], 6, 7)
  top.table <- as.data.frame(top.table)
  try(top.table$symbol <- mapIds(org.Hs.eg.db, keys = top.table$Geneid,
                                 keytype = "ENSEMBL", column = "SYMBOL",
                                 multiVals="first")) #adding gene names
  top.table <- subset(top.table, top.table$symbol != 'NA')
  top.table<- top.table %>%
    mutate(direction = case_when(logFC > 0.5 ~ "up",
                                 logFC < -0.5 ~ "down"))
  top.table <- top.table[(which(top.table$adj.P.Val <= 0.1 &
                                  abs(top.table$logFC) >= 0.5)),]
  d3Mb[[i]] = top.table$Geneid
  names(d3Mb)[i] <- exposure
  rownames(top.table) <- NULL
  top.table <- top.table [(c(7,1,2,3,4,5,6,8,9))]
  a3Mb[[i]] = top.table$symbol
  names(a3Mb)[i] <- exposure
  b3Mb[[i]] <- top.table
  names(b3Mb)[i] <- exposure}}


# removing the empty control group from the lists
a3Mb[3] <- NULL
b3Mb[3] <- NULL
d3Mb[3] <- NULL
```

``` r
###########################################################
p <- list()
p[[1]]<- a3Fb$CF
names(p)[1] <- "Female"

p[[2]]<- a3Mb$CF
names(p)[2] <- "Male"

EP_CF_FvM_3X <- euler(p, shape = "circle") # Euler plot 
plot(EP_CF_FvM_3X,
     quantities = TRUE, cex = 5,
     main = "Cardboard Flaming -repeated",
     lty = 1:3,
     labels = list(font = 3, cex = 2),
     fills = c( "khaki3", "white" ))# "burlywood1"
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-1.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fb$PF
names(p)[1] <- "Female"

p[[2]]<- a3Mb$PF
names(p)[2] <- "Male"

EP_CF_FvM_3X <- euler(p, shape = "circle") # Euler plot 
plot(EP_CF_FvM_3X,
     quantities = TRUE, cex = 5,
     main = "Plastic Flaming -repeated",
     lty = 1:3,
     labels = list(font = 3, cex = 2),
     fills = c( "khaki3", "white" ))
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-2.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fb$CS
names(p)[1] <- "Female"

p[[2]]<- a3Mb$CS
names(p)[2] <- "Male"


p <- ggVennDiagram(p, label = "count") + ggtitle("Cardboard Smoldering -repeated")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-3.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fb$PS
names(p)[1] <- "Female"

p[[2]]<- a3Mb$PS
names(p)[2] <- "Male"


p <- ggVennDiagram(p, label = "count") + ggtitle("Plastic Smoldering -repeated")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-4.png)<!-- -->

``` r
###########################################################
# ------------------------------------------------------ #
###########################################################
p <- list()
p[[1]]<- a3Fa$CF
names(p)[1] <- "Female"

p[[2]]<- a3Ma$CF
names(p)[2] <- "Male"

EP_CF_FvM_3X <- euler(p, shape = "circle") # Euler plot 
plot(EP_CF_FvM_3X,
     quantities = TRUE, cex = 5,
     main = "Cardboard Flaming -single",
     lty = 1:3,
     labels = list(font = 3, cex = 2),
     fills = c( "khaki3", "white" ))
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-5.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fa$PF
names(p)[1] <- "Female"

p[[2]]<- a3Ma$PF
names(p)[2] <- "Male"

EP_CF_FvM_3X <- euler(p, shape = "circle") # Euler plot 
plot(EP_CF_FvM_3X,
     quantities = TRUE, cex = 5,
     main = "Plastic Flaming -single",
     lty = 1:3,
     labels = list(font = 3, cex = 2),
     fills = c( "khaki3", "white" ))
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-6.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fa$CS
names(p)[1] <- "Female"

p[[2]]<- a3Ma$CS
names(p)[2] <- "Male"


p <- ggVennDiagram(p, label = "count") + ggtitle("Cardboard Smoldering -single")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-7.png)<!-- -->

``` r
###########################################################
p <- list()
p[[1]]<- a3Fa$PS
names(p)[1] <- "Female"

p[[2]]<- a3Ma$PS
names(p)[2] <- "Male"


p <- ggVennDiagram(p, label = "count") + ggtitle("Plastic Smoldering -single")+
  scale_x_continuous(expand = expansion(mult = .2))+
  ggplot2::scale_fill_gradient(low = "white", high = "white")
p
```

![](README_figs/README--Euler%20plots%20and%20Venn%20diagrams%20comapring%20Female%20vs%20Male-8.png)<!-- -->

``` r
##########################################################################

Gdata <- c3Fb # female repeated

#---------------------------------------------#
# Cardboard exposure group -Flaming

Gdata$CF <- Gdata$CF %>% mutate(ProbeID = Gdata$CF$Geneid)

GdataCBf <- Gdata$CF #Cardboard
GdataCBf <- select(GdataCBf, ProbeID, logFC)
rownames(GdataCBf) <- NULL

#making ranked gene list
genelist_GdataCBf = GdataCBf[,2] #numeric vector
names(genelist_GdataCBf) = as.character(GdataCBf[,1]) #named vector
genelist_GdataCBf = sort(genelist_GdataCBf, decreasing = TRUE) #must sort in descending order



#--------------------------------------------------------------------#
#--------------------------------------------------------------------#

#GO over-representation analysis
geneCBf <- names(genelist_GdataCBf[genelist_GdataCBf[] > 0.5]) # logFC > 0.5
enGOCBf <- enrichGO(gene         = geneCBf,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)

#head(enGOCBf)

enGOCBf <- gsfilter(enGOCBf, by = 'Count', min = 3) # selecting categories with at least gene count 3

#barplot
barplot(enGOCBf, showCategory=10, font = 15, title = "Cardboard Flaming -female")+
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))+
  scale_y_discrete(label=function(y) stringr::str_trunc(y, 40))
```

![](README_figs/README--GO%20enrichment%20of%20Female%20donors%20data-1.png)<!-- -->

``` r
##################################################################
#---------------------------------------------#
# Plastic exposure group -Flaming

Gdata$PF <- Gdata$PF %>% mutate(ProbeID = Gdata$PF$Geneid)   

GdataPLf <- Gdata$PF #Plastic
GdataPLf <- select(GdataPLf, ProbeID, logFC)
rownames(GdataPLf) <- NULL

#making ranked gene list
genelist_GdataPLf = GdataPLf[,2] #numeric vector
names(genelist_GdataPLf) = as.character(GdataPLf[,1]) #named vector
genelist_GdataPLf = sort(genelist_GdataPLf, decreasing = TRUE) #must sort in descending order

#--------------------------------------------------------------------#
#--------------------------------------------------------------------#
#--------------------------------------------------------------------#

#GO over-representation analysis
genePLf <- names(genelist_GdataPLf[genelist_GdataPLf[] > 0.5]) # logFC > 0.5
enGOPLf <- enrichGO(gene         = genePLf,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.1)

#head(enGOPLf)

enGOPLf <- gsfilter(enGOPLf, by = 'Count', min = 3) # selecting categories with at least gene count 3

#barplot
barplot(enGOPLf, showCategory=10, font = 15, title = "Plastic Flaming -female")+ 
  scale_x_continuous(breaks = seq(0, 20, by = 5), limits=c(0,20))+
  theme(axis.text.y = element_text(lineheight = 0.7, size = 15),
        title = element_text(size = 15, face="bold"))+
  scale_y_discrete(label=function(y) stringr::str_trunc(y, 40))
```

![](README_figs/README--GO%20enrichment%20of%20Female%20donors%20data-2.png)<!-- -->

``` r
#################################################################

k <- list()

k[[1]] <- d1$CF
names(k)[1] <- "sCF"

k[[2]] <- d1$PF
names(k)[2] <- "sPF"

k[[3]] <- d2$CF
names(k)[3] <- "mCF"

k[[4]] <- d2$PF
names(k)[4] <- "mPF"

z <- process_region_data(Venn(k))
z <- as.data.frame(z)

counts8 <- counts %>%
  select(contains("Geneid") | 
           contains("CF") | contains("PF") | contains("PBS") ) 

rownames(counts8) <- counts8$Geneid
counts8 <- counts8 %>% select(-c(Geneid) ) 
#-----------------------------------------------------------------------#
# Control and flaming
counts8 <- subset(counts8, rownames(counts8) %in% c(unlist(z[15,3])))

counts8$gene <- mapIds(org.Hs.eg.db, keys = row.names(counts8), 
                          keytype = "ENSEMBL", column = "SYMBOL", 
                          multiVals="first") #adding gene names 
counts8 <- subset(counts8, counts8$gene != 'NA')
rownames(counts8) <- NULL
rownames(counts8) <- counts8$gene
counts8 <- counts8 %>% select(-c(gene) ) 

x <- as.data.frame(colnames(counts8))
rownames(x) <- colnames(counts8)

colnames(x) <- "Groups"
x$Groups <- substr(x$Groups, 10, nchar(x$Groups))

x$Donors <- rownames(x)

x$Groups <- str_replace_all(x$Groups, c("PBS_Mock_1" = "Control_single",
                                        "PBS_Mock_3" = "Control_repeated",
                                        "PF_Mock_1" = "Plastic_single",
                                        "PF_Mock_3" = "Plastic_repeated",
                                        "CF_Mock_1" = "Cardboard_single",
                                        "CF_Mock_3" = "Cardboard_repeated"))

x$Donors <- substr(x$Donors, 6, 8)

x$Donors <- str_replace_all(x$Donors, c("15R" = "Female",
                                        "40R" = "Female",
                                        "68J" = "Male",
                                        "21O" = "Male",
                                        "67L" = "Male",
                                        "57N" = "Female"))

x2 = list( Groups = c(Cardboard_single = "green", Cardboard_repeated = "seagreen",
                      Plastic_single = "skyblue", Plastic_repeated = "navyblue",
                      Control_single = "lightyellow", Control_repeated = "yellow"),
           Donors = c(Female = "violet", 
                      Male = "purple"))

temp3 <- as.matrix(counts8)

pheatmap(temp3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), # sets color scheme
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5,
         cellwidth = 9, 
         cellheight = 9, 
         border_color = "black",
         annotation_col = x,
         annotation_colors = x2,
         show_colnames = F,
         main = 'genes altered by flaming condensates\nsingle and repeated exposures',
         treeheight_col = 9, 
         treeheight_row = 9,
         fontsize_row = 7, 
         scale = 'row', 
         fontsize_col = 7, 
         cutree_rows = 2, 
         cluster_cols = TRUE, 
         cluster_rows = TRUE)
```

![](README_figs/README--Heat%20map%20of%20common%20genes-1.png)<!-- -->

``` r
#------------------------------------------------------------#

# for common genes altered by smoke condensates -heat maps 
# (single and repeated exposures)

HM1 <- counts

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------#
HM1a <- HM1 %>% 
  select(contains("Geneid") | contains("_PBS_") | contains("_CF_"))

rownames(HM1a) <- HM1a$Geneid
HM1a <- HM1a %>% select(-c(Geneid) ) 

#-----------------------------------------------------------------------#

CBf <- list()

CBf[[1]]<- d1$CF
names(CBf)[1] <- "single"

CBf[[2]]<- d2$CF
names(CBf)[2] <- "multiple"

z <- process_region_data(Venn(CBf))
z <- as.data.frame(z)
z[3,3]

counts_ALL <- subset(HM1a, rownames(HM1a) %in% unlist(z[3,3])) #genes

counts_ALL$gene <- mapIds(org.Hs.eg.db, keys = row.names(counts_ALL), 
                          keytype = "ENSEMBL", column = "SYMBOL", 
                          multiVals="first") #adding gene names 
counts_ALL <- subset(counts_ALL, counts_ALL$gene != 'NA')
rownames(counts_ALL) <- NULL
rownames(counts_ALL) <- counts_ALL$gene
counts_ALL <- counts_ALL[,1:(ncol(counts_ALL)-1)]

x <- as.data.frame(colnames(counts_ALL))
rownames(x) <- colnames(counts_ALL)

colnames(x) <- "Groups"

x$Groups <- substr(x$Groups, 10, nchar(x$Groups))


colnames(x) <- "Groups"
x$Donors <- rownames(x)

x$Groups <- str_replace_all(x$Groups, c("PBS_Mock_1" = "Control_single",
                                        "PBS_Mock_3" = "Control_repeated",
                                        "CF_Mock_1" = "Flaming_single",
                                        "CF_Mock_3" = "Flaming_repeated"))

x$Donors <- substr(x$Donors, 6, 8)

x$Donors <- str_replace_all(x$Donors, c("15R" = "Female",
                                              "40R" = "Female",
                                              "68J" = "Male",
                                              "21O" = "Male",
                                              "67L" = "Male",
                                              "57N" = "Female"))

x2 = list( Groups = c(Flaming_repeated = "orange", Flaming_single = "yellow",
                      Control_repeated = "darkgray", Control_single = "white"),
           Donors = c(Female = "violet", 
                      Male = "purple"))


temp3 <- as.matrix(counts_ALL)

pheatmap(temp3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), # sets color scheme
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5, 
         cellwidth = 9, 
         cellheight = 7, 
         border_color = "black",
         annotation_col = x,
         annotation_colors = x2,
         show_colnames = F,
         main = 'Cardboard flaming\n(common between single and repeated exposures)',
         treeheight_col = 9, 
         fontsize_row = 7, 
         scale = 'row', 
         fontsize_col = 7, 
         cutree_rows = 2, 
         cluster_cols = TRUE, 
         cluster_rows = TRUE)
```

![](README_figs/README--heat%20maps-1.png)<!-- -->

``` r
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
HM1b <- HM1 %>% 
  select(contains("Geneid") | contains("_PBS_") | contains("_CS_"))

rownames(HM1b) <- HM1b$Geneid
HM1b <- HM1b %>% select(-c(Geneid) ) 

#-----------------------------------------------------------------------#

CBs <- list()

CBs[[1]]<- d1$CS
names(CBs)[1] <- "single"

CBs[[2]]<- d2$CS
names(CBs)[2] <- "multiple"

z <- process_region_data(Venn(CBs))
z <- as.data.frame(z)
z[3,3]

counts_ALL <- subset(HM1b, rownames(HM1b) %in% unlist(z[3,3])) #genes

counts_ALL$gene <- mapIds(org.Hs.eg.db, keys = row.names(counts_ALL), 
                          keytype = "ENSEMBL", column = "SYMBOL", 
                          multiVals="first") #adding gene names 
counts_ALL <- subset(counts_ALL, counts_ALL$gene != 'NA')
rownames(counts_ALL) <- NULL
rownames(counts_ALL) <- counts_ALL$gene
counts_ALL <- counts_ALL[,1:(ncol(counts_ALL)-1)]

x <- as.data.frame(colnames(counts_ALL))
rownames(x) <- colnames(counts_ALL)

colnames(x) <- "Groups"

x$Groups <- substr(x$Groups, 10, nchar(x$Groups))

colnames(x) <- "Groups"
x$Donors <- rownames(x)

x$Groups <- str_replace_all(x$Groups, c("PBS_Mock_1" = "Control_single",
                                        "PBS_Mock_3" = "Control_repeated",
                                        "CS_Mock_1" = "Smoldering_single",
                                        "CS_Mock_3" = "Smoldering_repeated"))

x$Donors <- substr(x$Donors, 6, 8)

x$Donors <- str_replace_all(x$Donors, c("15R" = "Female",
                                              "40R" = "Female",
                                              "68J" = "Male",
                                              "21O" = "Male",
                                              "67L" = "Male",
                                              "57N" = "Female"))

x2 = list( Groups = c(Smoldering_repeated = "orange", Smoldering_single = "yellow",
                      Control_repeated = "darkgray", Control_single = "white"),
           Donors = c(Female = "violet", 
                      Male = "purple"))


temp3 <- as.matrix(counts_ALL)

pheatmap(temp3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), # sets color scheme
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5, 
         cellwidth = 9, 
         cellheight = 7, 
         border_color = "black",
         annotation_col = x,
         annotation_colors = x2,
         show_colnames = F,
         main = 'Cardboard smoldering\n(common between single and repeated exposures)',
         treeheight_col = 9, 
         fontsize_row = 7, 
         scale = 'row', 
         fontsize_col = 7,
         cluster_cols = TRUE, 
         cluster_rows = TRUE)
```

![](README_figs/README--heat%20maps-2.png)<!-- -->

``` r
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
HM1c <- HM1 %>% 
  select(contains("Geneid") | contains("_PBS_") | contains("_PF_"))

rownames(HM1c) <- HM1c$Geneid
HM1c <- HM1c %>% select(-c(Geneid) ) 

#-----------------------------------------------------------------------#

PLf <- list()

PLf[[1]]<- d1$PF
names(PLf)[1] <- "single"

PLf[[2]]<- d2$PF
names(PLf)[2] <- "multiple"

z <- process_region_data(Venn(PLf))
z <- as.data.frame(z)
z[3,3]

counts_ALL <- subset(HM1c, rownames(HM1c) %in% unlist(z[3,3])) #genes

counts_ALL$gene <- mapIds(org.Hs.eg.db, keys = row.names(counts_ALL), 
                          keytype = "ENSEMBL", column = "SYMBOL", 
                          multiVals="first") #adding gene names 
counts_ALL <- subset(counts_ALL, counts_ALL$gene != 'NA')
rownames(counts_ALL) <- NULL
rownames(counts_ALL) <- counts_ALL$gene
counts_ALL <- counts_ALL[,1:(ncol(counts_ALL)-1)]

x <- as.data.frame(colnames(counts_ALL))
rownames(x) <- colnames(counts_ALL)

colnames(x) <- "Groups"

x$Groups <- substr(x$Groups, 10, nchar(x$Groups))

colnames(x) <- "Groups"
x$Donors <- rownames(x)

x$Groups <- str_replace_all(x$Groups, c("PBS_Mock_1" = "Control_single",
                                        "PBS_Mock_3" = "Control_repeated",
                                        "PF_Mock_1" = "Flaming_single",
                                        "PF_Mock_3" = "Flaming_repeated"))

x$Donors <- substr(x$Donors, 6, 8)

x$Donors <- str_replace_all(x$Donors, c("15R" = "Female",
                                        "40R" = "Female",
                                        "68J" = "Male",
                                        "21O" = "Male",
                                        "67L" = "Male",
                                        "57N" = "Female"))

x2 = list( Groups = c(Flaming_repeated = "orange", Flaming_single = "yellow",
                      Control_repeated = "darkgray", Control_single = "white"),
           Donors = c(Female = "violet", 
                      Male = "purple"))


temp3 <- as.matrix(counts_ALL)

pheatmap(temp3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), # sets color scheme
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5, 
         cellwidth = 9, 
         cellheight = 7, 
         border_color = "black",
         annotation_col = x,
         annotation_colors = x2,
         show_colnames = F,
         main = 'Plastic flaming\n(common between single and repeated exposures)',
         treeheight_col = 9, 
         fontsize_row = 7, 
         scale = 'row', 
         fontsize_col = 7, 
         cutree_rows = 2,
         cluster_cols = TRUE, 
         cluster_rows = TRUE)
```

![](README_figs/README--heat%20maps-3.png)<!-- -->

``` r
#-----------------------------------------------------------------------#
#-----------------------------------------------------------------------#
HM1d <- HM1 %>% 
  select(contains("Geneid") | contains("_PBS_") | contains("_PS_"))

rownames(HM1d) <- HM1d$Geneid
HM1d <- HM1d %>% select(-c(Geneid) ) 


#-----------------------------------------------------------------------#

PLs <- list()

PLs[[1]]<- d1$PS
names(PLs)[1] <- "single"

PLs[[2]]<- d2$PS
names(PLs)[2] <- "multiple"

z <- process_region_data(Venn(PLs))
z <- as.data.frame(z)
z[3,3]

counts_ALL <- subset(HM1d, rownames(HM1d) %in% unlist(z[3,3])) #genes

counts_ALL$gene <- mapIds(org.Hs.eg.db, keys = row.names(counts_ALL), 
                          keytype = "ENSEMBL", column = "SYMBOL", 
                          multiVals="first") #adding gene names 
counts_ALL <- subset(counts_ALL, counts_ALL$gene != 'NA')
rownames(counts_ALL) <- NULL
rownames(counts_ALL) <- counts_ALL$gene
counts_ALL <- counts_ALL[,1:(ncol(counts_ALL)-1)]

x <- as.data.frame(colnames(counts_ALL))
rownames(x) <- colnames(counts_ALL)

colnames(x) <- "Groups"

x$Groups <- substr(x$Groups, 10, nchar(x$Groups))

colnames(x) <- "Groups"
x$Donors <- rownames(x)

x$Groups <- str_replace_all(x$Groups, c("PBS_Mock_1" = "Control_single",
                                        "PBS_Mock_3" = "Control_repeated",
                                        "PS_Mock_1" = "Smoldering_single",
                                        "PS_Mock_3" = "Smoldering_repeated"))

x$Donors <- substr(x$Donors, 6, 8)

x$Donors <- str_replace_all(x$Donors, c("15R" = "Female",
                                        "40R" = "Female",
                                        "68J" = "Male",
                                        "21O" = "Male",
                                        "67L" = "Male",
                                        "57N" = "Female"))

x2 = list( Groups = c(Smoldering_repeated = "orange", Smoldering_single = "yellow",
                      Control_repeated = "darkgray", Control_single = "white"),
           Donors = c(Female = "violet", 
                      Male = "purple"))


temp3 <- as.matrix(counts_ALL)

pheatmap(temp3,
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), # sets color scheme
         display_numbers = FALSE,
         number_color = "black",
         fontsize_number = 5, 
         cellwidth = 9, 
         cellheight = 7, 
         border_color = "black",
         annotation_col = x,
         annotation_colors = x2,
         show_colnames = F,
         main = 'Plastic smoldering\n(common between single and repeated exposures)',
         treeheight_col = 9, 
         fontsize_row = 7, 
         scale = 'row', 
         fontsize_col = 7, 
         cluster_cols = TRUE, 
         cluster_rows = TRUE)
```

![](README_figs/README--heat%20maps-4.png)<!-- -->

``` r
#------------------------------------------------------------#
#-----------------------------------------------------------------------#
```

``` r
###### SAVING SUPPLEMENTARY TABLES WITH SIGNIFICANT DEGs ####
##############################################################

names(b1) <- str_replace_all(names(b1), c(CF= "Cardboard Flaming -single",
                                          CS= "Cardboard Smoldering -single",
                                          PF= "Plastic Flaming -single",
                                          PS= "Plastic Smoldering -single"))

names(b2) <- str_replace_all(names(b2), c(CF= "Cardboard Flaming -repeated",
                                          CS= "Cardboard Smoldering -repeated",
                                          PF= "Plastic Flaming -repeated",
                                          PS= "Plastic Smoldering -repeated"))

u1 <- c(b1, b2)

col_arrangement <- c("Geneid", "symbol",    "logFC", "AveExpr", "t", "B",   "P.Value", "adj.P.Val")

u1 <- map(u1, ~ relocate(.x, any_of(col_arrangement)))
u1 <- u1 %>% lapply(arrange, direction)
blank_excel <- createWorkbook()

Map(function(df, tab_name){     
  
  addWorksheet(blank_excel, tab_name)
  writeData(blank_excel, tab_name, df)
}, 

u1, names(u1)
)

saveWorkbook(blank_excel, file = "Supplementary Table 1.xlsx", overwrite = TRUE)

##############################################################

##############################################################
names(b3Fa) <- str_replace_all(names(b3Fa), c(CF= "Cardboard Flaming -single",
                                          CS= "Cardboard Smoldering -single",
                                          PF= "Plastic Flaming -single",
                                          PS= "Plastic Smoldering -single"))


names(b3Fb) <- str_replace_all(names(b3Fb), c(CF= "Cardboard Flaming -repeated",
                                              CS= "Cardboard Smoldering -repeated",
                                              PF= "Plastic Flaming -repeated",
                                              PS= "Plastic Smoldering -repeated"))

u2 <- c(b3Fa, b3Fb)

u2 <- map(u2, ~ relocate(.x, any_of(col_arrangement)))
u2 <- u2 %>% lapply(arrange, direction)
blank_excel <- createWorkbook()

Map(function(df, tab_name){     
  
  addWorksheet(blank_excel, tab_name)
  writeData(blank_excel, tab_name, df)
}, 

u2, names(u2)
)

saveWorkbook(blank_excel, file = "Supplementary Table 2.xlsx", overwrite = TRUE)

##############################################################
names(b3Ma) <- str_replace_all(names(b3Ma), c(CF= "Cardboard Flaming -single",
                                              CS= "Cardboard Smoldering -single",
                                              PF= "Plastic Flaming -single",
                                              PS= "Plastic Smoldering -single"))


names(b3Mb) <- str_replace_all(names(b3Mb), c(CF= "Cardboard Flaming -repeated",
                                              CS= "Cardboard Smoldering -repeated",
                                              PF= "Plastic Flaming -repeated",
                                              PS= "Plastic Smoldering -repeated"))

u3 <- c(b3Ma, b3Mb)

u3 <- map(u3, ~ relocate(.x, any_of(col_arrangement)))
u3 <- u3 %>% lapply(arrange, direction)
blank_excel <- createWorkbook()

Map(function(df, tab_name){     
  
  addWorksheet(blank_excel, tab_name)
  writeData(blank_excel, tab_name, df)
}, 

u3, names(u3)
)

saveWorkbook(blank_excel, file = "Supplementary Table 3.xlsx", overwrite = TRUE)
```

``` r
#getwd()
```
