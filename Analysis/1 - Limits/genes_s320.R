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
s320 <- c("ERR2098382","ERR2098383","ERR2098384", "ERR2098391", "ERR2098392","ERR2098393", "ERR2098402","ERR2098403","ERR2098404","ERR2098405", "ERR2098410","ERR2098411")
genes_s320 <- genes %>% select(c(1,all_of(s320)))
rm(genes)
genes_s320 <- t(genes_s320)

genes_s320_acc <- specaccum(genes_s320, method = "exact", permutations = 1000)

dev.new()
png("genes_s320_acc.png")
plot(genes_s320_acc)
dev.off()