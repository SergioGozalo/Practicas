library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

onel <- c("ERR2098365","ERR2098366","ERR2098367")
genes_onel <- genes %>% select(c(0,all_of(onel)))
genes_onel <- t(genes_onel)

genes_onel_acc <- specaccum(genes_onel, method = "exact", permutations = 100)

dev.new()
png("genes_onel_acc.png")
plot(genes_onel_acc)
dev.off()