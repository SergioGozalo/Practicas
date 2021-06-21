library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

s02 <- c("ERR2098375","ERR2098376","ERR2098377","ERR2098365","ERR2098366","ERR2098367","ERR2098368","ERR2098369","ERR2098370","ERR2098371","ERR2098372","ERR2098373","ERR2098374")

genes_s02_s <- genes %>% select(c(0,all_of(s02)))
genes_s02_s <- t(genes_s02_s)

genes_s02_s_acc <- specaccum(genes_s02_s, method = "exact", permutations = 100)

specslope(genes_s02_s_acc, 12 )
dev.new()
png("genes_s02_s_acc.png")
plot(genes_s02_s_acc, main = ">0.2")
dev.off()
