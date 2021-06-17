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
s23 <- c("ERR2098379","ERR2098380","ERR2098381", "ERR2098388","ERR2098389","ERR2098390", "ERR2098398","ERR2098399","ERR2098400", "ERR2098401", "ERR2098408", "ERR2098409", "ERR2196984")
genes_s23 <- genes %>% select(c(1,all_of(s23)))
rm(genes)
genes_s23 <- t(genes_s23)

genes_s23_acc <- specaccum(genes_s23, method = "exact", permutations = 1000)

dev.new()
png("genes_s23_acc.png")
plot(genes_s23_acc)
dev.off()