setwd("/Users/Usuario/Desktop/Bioinformatica/Practicas/Practicas/Analysis/"))


## Loading necessary libraries


library("vegan")
library("tidyr")
library("plyr")
library("dplyr")
library("reshape2")
library("tidyverse")


## Table reading


# Solo un fragmento de la tabla para ver que funciona todo
genes <- read.delim("functional.tables/prueba3.tbl.txt", header = TRUE, sep = "\t", row.names = NULL)


## Genetic diversity


gennames <- genes$gen
genes$gene<- NULL

# Hay que cambiar las comas de los exponentes por puntos, sino da error
genes <- mutate_if(genes, 
                   is.character, 
                   str_replace_all, pattern = ",", replacement = ".")

# Transponer, transformar en df y en num
genes <- as.data.frame(sapply(genes, as.numeric))
genes2 <- genes
rownames(genes2) <- gennames
genes2 <- t(genes2)



#Rarecurve
gene_acc <- specaccum(genes2, method = "exact", permutations = 1000)

dev.new()
png("genetic_diversity.png")
plot(gene_acc)
dev.off()



## By filter


# Select samples by filter
s02_s <- c("ERR2098365","ERR2098366","ERR2098367","ERR2098368","ERR2098369","ERR2098370","ERR2098371","ERR2098372","ERR2098373","ERR2098374")
s02_m <- c("ERR2098375","ERR2098376","ERR2098377")
s20 <- c("ERR2098378", "ERR2098385","ERR2098386","ERR2098387", "ERR2098394", "ERR2098395","ERR2098396","ERR2098397", "ERR2098406","ERR2098407", "ERR2196983")
s23 <- c("ERR2098379","ERR2098380","ERR2098381", "ERR2098388","ERR2098389","ERR2098390", "ERR2098398","ERR2098399","ERR2098400", "ERR2098401", "ERR2098408", "ERR2098409", "ERR2196984")
s320 <- c("ERR2098382","ERR2098383","ERR2098384", "ERR2098391", "ERR2098392","ERR2098393", "ERR2098402","ERR2098403","ERR2098404","ERR2098405", "ERR2098410","ERR2098411")

# Create the dataset for each filter and traspose
genes_s02_s <- genes %>% select(c(1,all_of(s02_s)))
genes_s02_s <- t(genes_s02_s)
genes_s02_m <- genes %>% select(c(1,all_of(s02_m)))
genes_s02_m <- t(genes_s02_m)
genes_s20 <- genes %>% select(c(1,all_of(s20)))
genes_s20 <- t(genes_s20)
genes_s23 <- genes %>% select(c(1,all_of(s23)))
genes_s23 <- t(genes_s23)
genes_s320 <- genes %>% select(c(1,all_of(s320)))
genes_s320 <- t(genes_s320)

# Accumulation curves
genes_s02_s_acc <- specaccum(genes_s02_s, method = "exact", permutations = 1000)

genes_s02_m_acc <- specaccum(genes_s02_m, method = "exact", permutations = 1000)

genes_s20_acc <- specaccum(genes_s20, method = "exact", permutations = 1000)

genes_s23_acc <- specaccum(genes_s23, method = "exact", permutations = 1000)

genes_s320_acc <- specaccum(genes_s320, method = "exact", permutations = 1000)

#Saving plots
dev.new()
png("genes_s02_s_acc.png")
plot(genes_s02_s_acc)
dev.off()

dev.new()
png("genes_s02_m_acc.png")
plot(genes_s02_m_acc)
dev.off()

dev.new()
png("genes_s20_acc.png")
plot(genes_s20_acc)
dev.off()

dev.new()
png("genes_s23_acc.png")
plot(genes_s23_acc)
dev.off()

dev.new()
png("genes_s23_acc.png")
plot(genes_s320_acc)
dev.off()

## By water volume


# Select samples by water volume
onel <- c("ERR2098365","ERR2098366","ERR2098367")
ten <- c("ERR2098373","ERR2098374","ERR2098375","ERR2098376","ERR2098377","ERR2098378","ERR2098379", "ERR2098380","ERR2098381","ERR2098382", "ERR2098383","ERR2098384", "ERR2098369","ERR2098370","ERR2098371","ERR2098372")
hundred <- c("ERR2098385","ERR2098386","ERR2098387","ERR2098388","ERR2098389","ERR2098390", "ERR2098391","ERR2098392", "ERR2098393", "ERR2098395","ERR2098396","ERR2098397", "ERR2098399","ERR2098400","ERR2098401", "ERR2098403","ERR2098404","ERR2098405", "ERR2098407","ERR2098409","ERR2098411")
five_hundred <- c("ERR2098394", "ERR2098398", "ERR2098402")
thousand <- c( "ERR2098406", "ERR2098408", "ERR2098410", "ERR2196983","ERR2196984", "ERR2196985")

# Create the dataset for each filter and traspose
genes_onel <- genes %>% select(c(1,all_of(onel)))
genes_onel <- t(genes_onel)
genes_ten <- genes %>% select(c(1,all_of(ten)))
genes_ten <- t(genes_ten)
genes_hundred <- genes %>% select(c(1,all_of(hundred)))
genes_hundred <- t(genes_hundred)
genes_five_hundred <- genes %>% select(c(1,all_of(five_hundred)))
genes_five_hundred <- t(genes_five_hundred)
genes_thousand <- genes %>% select(c(1,all_of(thousand)))
genes_thousand <- t(genes_thousand)


# Accumulation curves
genes_onel_acc <- specaccum(genes_onel, method = "exact", permutations = 1000)

genes_five_hundred_acc <- specaccum(genes_five_hundred, method = "exact", permutations = 1000)

genes_ten_acc <- specaccum(genes_ten, method = "exact", permutations = 1000)

genes_hundred_acc <- specaccum(genes_hundred, method = "exact", permutations = 1000)

genes_thousand_acc <- specaccum(genes_thousand, method = "exact", permutations = 1000)



#Saving plots
dev.new()
png("genes_onel_acc.png")
plot(genes_onel_acc)
dev.off()

dev.new()
png("genes_five_hundred_acc.png")
plot(genes_five_hundred_acc)
dev.off()

dev.new()
png("genes_ten_acc.png")
plot(genes_ten_acc)
dev.off()

dev.new()
png("genes_hundred_acc.png")
plot(genes_hundred_acc)
dev.off()

dev.new()
png("genes_thousand_acc.png")
plot(genes_thousand_acc)
dev.off()
