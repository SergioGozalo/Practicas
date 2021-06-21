library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

#ten <- c("ERR2098365","ERR2098366","ERR2098367","ERR2098373","ERR2098374","ERR2098375","ERR2098376","ERR2098377","ERR2098378","ERR2098379", "ERR2098380","ERR2098381","ERR2098382", "ERR2098383","ERR2098384", "ERR2098369","ERR2098370","ERR2098371","ERR2098372")
ten <- c("ERR2098365","ERR2098366","ERR2098367","ERR2098373","ERR2098374","ERR2098375","ERR2098376","ERR2098377","ERR2098378","ERR2098379")

genes_ten <- genes %>% select(c(0,all_of(ten)))

genes_ten <- t(genes_ten)

genes_ten_acc <- specaccum(genes_ten, method = "exact", permutations = 100)

dev.new()
png("genes_ten_acc2.png")
plot(genes_ten_acc, main = "Ten liters")
dev.off()