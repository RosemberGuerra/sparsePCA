### Title:     Big5 data set.
### Created:  22-01-2020
### Modified: 

library(qgraph)


current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()

rm(list = ls(all.names = TRUE))

# loading data set #
data("big5")

write.table(big5,file = "big5.txt")
