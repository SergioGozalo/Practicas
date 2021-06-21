library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

s23 <- c("ERR2098379","ERR2098380","ERR2098381", "ERR2098388","ERR2098389","ERR2098390", "ERR2098398","ERR2098399","ERR2098400", "ERR2098401", "ERR2098408", "ERR2098409", "ERR2196984")
genes_s23 <- genes %>% select(c(0,all_of(s23)))
genes_s23 <- t(genes_s23)

genes_s23_acc <- specaccum(genes_s23, method = "exact", permutations = 100)

specslope(genes_s23_acc, 12 )
dev.new()
png("genes_s23_acc.png")
plot(genes_s23_acc, main = "0.22-3")
dev.off()