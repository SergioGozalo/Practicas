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
s02_s <- c("ERR2098365","ERR2098366","ERR2098367","ERR2098368","ERR2098369","ERR2098370","ERR2098371","ERR2098372","ERR2098373","ERR2098374")

genes_s02_s <- genes %>% select(c(1,all_of(s02_s)))
rm(genes)
genes_s02_s <- t(genes_s02_s)

genes_s02_s_acc <- specaccum(genes_s02_s, method = "exact", permutations = 1000)

dev.new()
png("genes_s02_s_acc.png")
plot(genes_s02_s_acc)
dev.off()
