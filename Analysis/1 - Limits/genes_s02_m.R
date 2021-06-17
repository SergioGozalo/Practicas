library("vegan")
library("tidyr")
library("plyr")
library("dplyr")
library("reshape2")
library("tidyverse")


genes <- read.delim("functional.tables/prueba3.tbl.txt", header = TRUE, sep = "\t", row.names = NULL)

gennames <- genes$gen
genes$gene<- NULL

# Hay que cambiar las comas de los exponentes por puntos, sino da error
genes <- mutate_if(genes, 
                   is.character, 
                   str_replace_all, pattern = ",", replacement = ".")

# Transponer, transformar en df y en num
genes <- as.data.frame(sapply(genes, as.numeric))
s02_m <- c("ERR2098375","ERR2098376","ERR2098377")
genes_s02_m <- genes %>% select(c(1,all_of(s02_m)))
rm(genes)
genes_s02_m <- t(genes_s02_m)

genes_s02_m_acc <- specaccum(genes_s02_m, method = "exact", permutations = 1000)

dev.new()
png("genes_s02_m_acc.png")
plot(genes_s02_m_acc)
dev.off()
