## Index of sparseness vs Proportion of sparsity #
## Gene expression ADNI data set.
## author: Rosember Guerra-Urzola
## created: 27-01-2020
## modified: 


library(ggplot2) # must be version > 2.2.0
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)
getwd()
rm(list = ls())
Gpower = read.csv("Result-Gpower_3.txt", sep = ",", header = TRUE)

Width = 12.50; Height = 8.94 #
# Varimax plot #
p <- ggplot(Gpower, aes(x = PS_Gpower))
p <- p + geom_line(aes(y = IS_Gpower, linetype = "IS"))
p <- p + geom_line(aes(y = PEV_Gpower/3, linetype = "PEV"))
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*3, name = "PEV [%]"))
p <- p + scale_linetype_manual(values=c("solid", "dashed"))
# p <- p + scale_fill_manual(values=c("#484848", "#808080")) +
# p <- p + geom_vline(xintercept = 0.9514, color = "red")
# p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Index of sparseness",
              x = "Proportion of Sparsity", linetype = " ",title = "Gpower")#+ ggtitle("Simplimax") 
# p <- p + theme(legend.position = c(0.8, 0.9))
p +  theme_bw()+ theme(plot.title = element_text(size=20 ,hjust = 0.5),
                       legend.position = c(0.2, 0.8),
                       axis.text=element_text(size=16),
                       axis.title = element_text(size = 20))
ggsave("IS_GPower.eps",width = Width, height = Height,
       limitsize = FALSE)
