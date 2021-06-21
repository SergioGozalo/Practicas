library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

five_hundred <- c("ERR2098394", "ERR2098398", "ERR2098402")
genes_five_hundred <- genes %>% select(c(0,all_of(five_hundred)))

genes_five_hundred <- t(genes_five_hundred)

genes_five_hundred_acc <- specaccum(genes_five_hundred, method = "exact", permutations = 100)
dev.new()
png("genes_five_hundred_acc.png")
plot(genes_five_hundred_acc)
dev.off()