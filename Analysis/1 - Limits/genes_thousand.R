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
thousand <- c( "ERR2098406", "ERR2098408", "ERR2098410", "ERR2196983","ERR2196984", "ERR2196985")
genes_thousand <- genes %>% select(c(1,all_of(thousand)))
rm(genes)
genes_thousand <- t(genes_thousand)

genes_thousand_acc <- specaccum(genes_thousand, method = "exact", permutations = 1000)
dev.new()
png("genes_thousand_acc.png")
plot(genes_thousand_acc)
dev.off()