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
ten <- c("ERR2098373","ERR2098374","ERR2098375","ERR2098376","ERR2098377","ERR2098378","ERR2098379", "ERR2098380","ERR2098381","ERR2098382", "ERR2098383","ERR2098384", "ERR2098369","ERR2098370","ERR2098371","ERR2098372")
genes_ten <- genes %>% select(c(1,all_of(ten)))
rm(genes)
genes_ten <- t(genes_ten)

genes_ten_acc <- specaccum(genes_ten, method = "exact", permutations = 1000)

dev.new()
png("genes_ten_acc.png")
plot(genes_ten_acc)
dev.off()