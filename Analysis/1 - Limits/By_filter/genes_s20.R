library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)


s20 <- c("ERR2098378", "ERR2098385","ERR2098386","ERR2098387", "ERR2098394", "ERR2098395","ERR2098396","ERR2098397", "ERR2098406","ERR2098407", "ERR2196983")
genes_s20 <- genes %>% select(c(0,all_of(s20)))

genes_s20 <- t(genes_s20)

genes_s20_acc <- specaccum(genes_s20, method = "exact", permutations = 100)

dev.new()
png("genes_s20_acc.png")
plot(genes_s20_acc, main = ">20")
dev.off()