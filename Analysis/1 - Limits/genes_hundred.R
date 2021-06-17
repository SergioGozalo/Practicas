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
hundred <- c("ERR2098385","ERR2098386","ERR2098387","ERR2098388","ERR2098389","ERR2098390", "ERR2098391","ERR2098392", "ERR2098393", "ERR2098395","ERR2098396","ERR2098397", "ERR2098399","ERR2098400","ERR2098401", "ERR2098403","ERR2098404","ERR2098405", "ERR2098407","ERR2098409","ERR2098411")
genes_hundred <- genes %>% select(c(1,all_of(hundred)))
rm(genes)
genes_hundred <- t(genes_hundred)

genes_hundred_acc <- specaccum(genes_hundred, method = "exact", permutations = 1000)
dev.new()
png("genes_hundred_acc.png")
plot(genes_hundred_acc)
dev.off()