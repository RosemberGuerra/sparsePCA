# Biplots and IS for Big 5 application #
# Rosember Guerra #
# 17-02-2020 #

library(R.matlab)
library(qgraph)
library(ggfortify)
library(ggplot2) # must be version > 2.2.0


# setting working directory #
current_working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(current_working_dir)

rm(list = ls(all.names = TRUE))

# loading the data sets #
data("big5")
sPCArSVD = readMat("sPCArSVDApplication.mat")

# class(sPCArSVD$X)
# View(sPCArSVD$X)
# View(big5)


# Plot- ordinary PCA # 
# dim ->  width=1071, height = 858 #

X = scale(big5)
pcX = svd(X)

PCscores = pcX$u[,1:5]
PCloadings = pcX$v[,1:5]%*%diag(pcX$d[1:5])

par(mfrow=c(1,1),pty="s")
plot(PCscores[,c(2,1)],pch=16,cex=0.5, cex.lab =2, ylab = "C. scores 1", xlab = "C. scores 2")
par(new=T,pty="s")
plot(PCloadings[,c(2,1)],col=4,axes=FALSE, pch="",xlab = " ",ylab = " ")
arrows (rep(0,nrow(PCloadings)),rep(0,nrow(PCloadings)),
        PCloadings[,2],PCloadings[,1],col="#808080",length=0.1,lwd=1)

axis (3)
axis (4) 

### plot- sparse PCA: sPCA-rSVD ###
# dim ->  width=1071, height = 858 #

indx = c(apply(sPCArSVD$P[,c(1,2)], 2, which.max),apply(sPCArSVD$P[,c(1,2)], 2, which.min))

par(mfrow=c(1,1),pty="s")
plot(sPCArSVD$T[,c(2,1)],pch=16,cex=0.5,cex.lab =2, ylab = "C. scores 1", xlab = "C. scores 2")#, col = Color2)

par(new=T,pty="s")
Rot <- sPCArSVD$P[,c(2,1)]

plot(Rot,col=4,axes=FALSE,#xlim=c(XLIM[1]-0.3,XLIM[2]),ylim=c(XLIM[1]-0.3,XLIM[2]),
     pch="",xlab = " ",ylab = " ")
text(x=Rot[indx,1],y=Rot[indx,2], colnames(big5)[indx], cex=0.65, pos=1)
arrows (rep(0,nrow(sPCArSVD$P)),rep(0,nrow(sPCArSVD$P)),Rot[,1],Rot[,2],col="#808080",length=0.1,lwd=1)

axis (3)
axis (4) 

###--- IS --- ###
sPCArSVD2 = read.csv("Result-sPCArSVD.txt", sep = ",", header = TRUE)

# sPCA-rSVD plot #
p <- ggplot(sPCArSVD2, aes(x = PS))
p <- p + geom_line(aes(y = IS_spcarsvd, linetype = "IS"))
p <- p + geom_line(aes(y = PEV_spcarsvd/5, linetype = "PEV"))
p <- p + scale_y_continuous(sec.axis = sec_axis(~.*5, name = "PEV"))
p <- p + scale_linetype_manual(values=c("solid", "dashed"))
# p <- p + scale_colour_manual(values = c("blue", "red"))
p <- p + labs(y = "Index of sparseness",
              x = "Proportion of Sparsity", linetype = " ",title = "sPCA-rSVD")#+ ggtitle("Simplimax") 
# p <- p + theme(legend.position = c(0.8, 0.9))
p +  theme_bw()+ theme(plot.title = element_text(size=20 ,hjust = 0.5),
                       legend.position = c(0.2, 0.8),
                       axis.text=element_text(size=16),
                       axis.title = element_text(size = 20))

# ggsave("sPCArSVD.eps",width = Width, height = Height,
#        limitsize = FALSE) 