library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

s320 <- c("ERR2098382","ERR2098383","ERR2098384", "ERR2098391", "ERR2098392","ERR2098393", "ERR2098402","ERR2098403","ERR2098404","ERR2098405", "ERR2098410","ERR2098411")
genes_s320 <- genes %>% select(c(1,all_of(s320)))
genes_s320 <- t(genes_s320)

genes_s320_acc <- specaccum(genes_s320, method = "exact", permutations = 100)

dev.new()
png("genes_s320_acc.png")
plot(genes_s320_acc, main ="3-20")
dev.off()