library("vegan")
library("tidyverse")


genes <- read.delim("functional.tables/table.txt", header = TRUE, sep = "\t", row.names = NULL)

#hundred <- c("ERR2098385","ERR2098386","ERR2098387","ERR2098388","ERR2098389","ERR2098390", "ERR2098391","ERR2098392", "ERR2098393", "ERR2098395","ERR2098396","ERR2098397", "ERR2098399","ERR2098400","ERR2098401", "ERR2098403","ERR2098404","ERR2098405", "ERR2098407","ERR2098409","ERR2098411")
hundred <- c("ERR2098385","ERR2098386","ERR2098387","ERR2098388","ERR2098389","ERR2098390", "ERR2098391","ERR2098392", "ERR2098393", "ERR2098395")

genes_hundred <- genes %>% select(c(0,all_of(hundred)))

genes_hundred <- t(genes_hundred)

genes_hundred_acc <- specaccum(genes_hundred, method = "exact", permutations = 100)
dev.new()
png("genes_hundred_acc2.png")
plot(genes_hundred_acc, main = "Hundred liters")
dev.off()