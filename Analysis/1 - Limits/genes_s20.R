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
s20 <- c("ERR2098378", "ERR2098385","ERR2098386","ERR2098387", "ERR2098394", "ERR2098395","ERR2098396","ERR2098397", "ERR2098406","ERR2098407", "ERR2196983")
genes_s20 <- genes %>% select(c(1,all_of(s20)))
rm(genes)
genes_s20 <- t(genes_s20)

genes_s20_acc <- specaccum(genes_s20, method = "exact", permutations = 1000)

dev.new()
png("genes_s20_acc.png")
plot(genes_s20_acc)
dev.off()