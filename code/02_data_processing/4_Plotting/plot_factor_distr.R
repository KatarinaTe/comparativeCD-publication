##################################################
####plot distribution using the factor scores ####
##################################################
famQC6_F1 <- read.csv("CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
famQC6_F2 <- read.csv("CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))
famQC6_F3 <- read.csv("CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","pheno"))

#plot distribution
(summary(famQC6_F1$pheno)[6])/(sd(famQC6_F1$pheno))
sd(famQC6_F1$pheno)
SD2VÄRDE <- sd(famQC6_F1$pheno)*2+summary(famQC6_F1$pheno)[4]
Q3VÄRDE <- summary(famQC6_F1$pheno)[5]
MEAN <- summary(famQC6_F1$pheno)[4]
MEDIAN <- summary(famQC6_F1$pheno)[3]
plotF1 <- ggplot(famQC6_F1, aes(x=pheno)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  annotate("text", x = Inf, y = Inf, label = "N=2429", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCDF1", x = "F1 score", y = "Number of dogs")+
  xlim(-2,3.5)+
  ylim(0,800) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + 
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        #axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              #fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))



(summary(famQC6_F2$pheno)[6])/(sd(famQC6_F2$pheno))
sd(famQC6_F2$pheno)
SD2VÄRDE <- sd(famQC6_F2$pheno)*2+summary(famQC6_F2$pheno)[4]
Q3VÄRDE <- summary(famQC6_F2$pheno)[5]
MEAN <- summary(famQC6_F2$pheno)[4]
MEDIAN <- summary(famQC6_F2$pheno)[3]
plotF2 <- ggplot(famQC6_F2, aes(x=pheno)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  annotate("text", x = Inf, y = Inf, label = "N=2506", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCDF2", x = "F2 score", y = "Number of dogs")+
  xlim(-2,3.5)+
  ylim(0,800) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + 
  theme(plot.title.position = "plot", 
      plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                          fill = "cornsilk"),
      #axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                            #fill = "white"),
      axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                            padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))


(summary(famQC6_F3$pheno)[6])/(sd(famQC6_F3$pheno))
sd(famQC6_F3$pheno)
SD2VÄRDE <- sd(famQC6_F3$pheno)*2+summary(famQC6_F3$pheno)[4]
Q3VÄRDE <- summary(famQC6_F3$pheno)[5]
MEAN <- summary(famQC6_F3$pheno)[4]
MEDIAN <- summary(famQC6_F3$pheno)[3]
plotF3 <- ggplot(famQC6_F3, aes(x=pheno)) + geom_histogram(color = "black", fill = "green", alpha = 0.5, binwidth = 0.1) +
  annotate("text", x = Inf, y = Inf, label = "N=2428", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCDF3", x = "F3 score", y = "Number of dogs")+
  xlim(-2,3.5)+
  ylim(0,800) +
  geom_vline(aes(xintercept = SD2VÄRDE), color = "red", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEAN), color = "purple", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = MEDIAN), color = "orange", linetype = "dashed", size = 0.5)+
  geom_vline(aes(xintercept = Q3VÄRDE), color = "blue", linetype = "dashed", size = 0.5) + 
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))


plotF1/plotF2/plotF3