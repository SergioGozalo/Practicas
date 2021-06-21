library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

thousand <- c( "ERR2098406", "ERR2098408", "ERR2098410", "ERR2196983","ERR2196984", "ERR2196985", "ERR2098394", "ERR2098398", "ERR2098402")
genes_thousand <- genes %>% select(c(0,all_of(thousand)))

genes_thousand <- t(genes_thousand)

genes_thousand_acc <- specaccum(genes_thousand, method = "exact", permutations = 100)
dev.new()
png("genes_thousand_acc.png")
plot(genes_thousand_acc, main = "Thousand liters")
dev.off()