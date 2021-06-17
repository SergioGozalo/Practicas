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
onel <- c("ERR2098365","ERR2098366","ERR2098367")
genes_onel <- genes %>% select(c(1,all_of(onel)))
rm(genes)
genes_onel <- t(genes_onel)

genes_onel_acc <- specaccum(genes_onel, method = "exact", permutations = 1000)

dev.new()
png("genes_onel_acc.png")
plot(genes_onel_acc)
dev.off()