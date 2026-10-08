#merge fam files:
CCDF1_fam <- read.csv("CCDF1_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","F1"))
CCDF2_fam <- read.csv("CCDF2_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","F2"))
CCDF3_fam <- read.csv("CCDF3_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","F3"))

CCDitem7_fam <- read.csv("CCDitem7_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item7"))
CCDitem93_fam <- read.csv("CCDitem93_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item93"))
CCDitem95_fam <- read.csv("CCDitem95_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item95"))
CCDitem145_fam <- read.csv("CCDitem145_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item145"))
CCDitem146_fam <- read.csv("CCDitem146_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item146"))
CCDitem147_fam <- read.csv("CCDitem147_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item147"))
CCDitem148_fam <- read.csv("CCDitem148_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item148"))
CCDitem149_fam <- read.csv("CCDitem149_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item149"))
CCDitem150_fam <- read.csv("CCDitem150_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item150"))
CCDitem151_fam <- read.csv("CCDitem151_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item151"))
CCDitem152_fam <- read.csv("CCDitem152_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item152"))
CCDitem153_fam <- read.csv("CCDitem153_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item153"))
CCDitem154_fam <- read.csv("CCDitem154_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item154"))
CCDitem155_fam <- read.csv("CCDitem155_DA_MERGED_GENCOVE_AXIOM_QC6.fam", sep= " ", header=F, col.names=c("FID","IID","F","M","sex","item155"))


F1.hist <- ggplot(CCDF1_fam, aes(x=F1)) +
  geom_histogram(alpha = 2,binwidth=0.1, color="darkseagreen3", fill="darkseagreen3") +
  scale_fill_manual(values=c("darkseagreen3")) +
  scale_color_manual(values=c("darkseagreen3")) +
  ylim(0,410) +
  scale_x_continuous(breaks = -1.5:3.5, limits = c(-1.6, 3.5)) +
  annotate("text", x = Inf, y = Inf, label = "N=2429", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCD factor 1: Repetition severity", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2,15), fill = "white"))


F2.hist <- ggplot(CCDF2_fam, aes(x=F2)) +
  geom_histogram(alpha = 2,binwidth=0.1, color="coral1", fill="coral1") +
  scale_fill_manual(values=c("coral1")) +
  scale_color_manual(values=c("coral1")) +
  ylim(0,410) +
  scale_x_continuous(breaks = -1.5:3.5, limits = c(-1.6, 3.5)) +
  annotate("text", x = Inf, y = Inf, label = "N=2506", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCD factor 2: Compulsive staring/trancing/pacing", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

F3.hist <- ggplot(CCDF3_fam, aes(x=F3)) +
  geom_histogram(alpha = 2,binwidth=0.1, color="dodgerblue2", fill="dodgerblue2") +
  scale_fill_manual(values=c("dodgerblue2")) +
  scale_color_manual(values=c("dodgerblue2")) +
  ylim(0,410) +
  scale_x_continuous(breaks = -1.5:3.5, limits = c(-1.6, 3.5)) +
  annotate("text", x = Inf, y = Inf, label = "N=2428", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "CCD factor 3: Repetition frequency", x = "Factor scores", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 38), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 0.65, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 15), fill = "white"))





F1.hist+F2.hist+F3.hist
summary(CCDF1_fam$F1, exclude=NULL)
summary(CCDF2_fam$F2, exclude=NULL)
summary(CCDF3_fam$F3, exclude=NULL)



####hist items for CCDF1####
library(patchwork)
library(ggtext)

table(CCDitem7_fam$item7)

item7.hist1 <- ggplot(CCDitem7_fam, aes(x=item7)) +
  geom_histogram(alpha = 1,binwidth=1, color="darkseagreen3", fill="darkseagreen3") +
  scale_fill_manual(values=c("darkseagreen3")) +
  scale_color_manual(values=c("darkseagreen3")) +
  ylim(0,2300) +
  annotate("text", x = Inf, y = Inf, label = "N=2555", hjust = 1.3, vjust = 2.5, size = 4) +
  labs(title = "item7: Excitement can lead DOG to fixed repetitive behavior", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 20), fill = "white"))


item153.hist1 <- ggplot(CCDitem153_fam, aes(x=item153)) +
  geom_histogram(alpha = 1,binwidth=1, color="darkslateblue", fill="darkslateblue") +
  scale_fill_manual(values=c("darkslateblue")) +
  scale_color_manual(values=c("darkslateblue")) +
  annotate("text", x = Inf, y = Inf, label = "N=2417", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item153: It can be difficult to get DOG's attention when HE is engaged in a repetitive behavior", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 0), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

item154.hist1 <- ggplot(CCDitem154_fam, aes(x=item154)) +
  geom_histogram(alpha = 1,binwidth=1, color="coral1", fill="coral1") +
  scale_fill_manual(values=c("coral1")) +
  scale_color_manual(values=c("coral1")) +
  annotate("text", x = Inf, y = Inf, label = "N=2418", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item154: It can be difficult to interrupt DOG when HE is engaged in a repetitive behavior", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 20), fill = "white"))


item155.hist1 <- ggplot(CCDitem155_fam, aes(x=item155)) +
  geom_histogram(alpha = 1,binwidth=1, color="dodgerblue2", fill="dodgerblue2") +
  scale_fill_manual(values=c("dodgerblue2")) +
  scale_color_manual(values=c("dodgerblue2")) +
  annotate("text", x = Inf, y = Inf, label = "N=2322", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item155: One or more repetitive behaviors interfere with DOG's life", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 0), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())




item7.hist1+item153.hist1+item154.hist1+item155.hist1


####hist items for CCDF2####

item93.hist1 <- ggplot(CCDitem93_fam, aes(x=item93)) +
  geom_histogram(alpha = 1,binwidth=1, color="turquoise3", fill="turquoise3") +
  scale_fill_manual(values=c("turquoise3")) +
  scale_color_manual(values=c("turquoise3")) +
  annotate("text", x = Inf, y = Inf, label = "N=2494", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item93: DOG paces up and down, walks in circles and/or wanders with no direction or purpose", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1.5, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))


item95.hist1 <- ggplot(CCDitem95_fam, aes(x=item95)) +
  geom_histogram(alpha = 1,binwidth=1, color="thistle3", fill="thistle3") +
  scale_fill_manual(values=c("thistle3")) +
  scale_color_manual(values=c("thistle3")) +
  annotate("text", x = Inf, y = Inf, label = "N=2498", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item95: DOG stares blankly at the walls or floor", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())


item150.hist1 <- ggplot(CCDitem150_fam, aes(x=item150)) +
  geom_histogram(alpha = 1,binwidth=1, color="dodgerblue4", fill="dodgerblue4") +
  scale_fill_manual(values=c("dodgerblue4")) +
  scale_color_manual(values=c("dodgerblue4")) +
  scale_x_continuous(breaks = 1:6) +
  annotate("text", x = Inf, y = Inf, label = "N=2417", hjust = 1.3, vjust = 2.5, size = 4) +
  ylim(0,2300) +
  labs(title = "item150: How much time in a normal day does DOG spend showing trance-like behaviour?", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())




item93.hist1+item95.hist1+item150.hist1

####hist items for CCDF3####

item145.hist1 <- ggplot(CCDitem145_fam, aes(x=item145)) +
  geom_histogram(alpha = 1,binwidth=1, color="turquoise3", fill="turquoise3") +
  scale_fill_manual(values=c("turquoise3")) +
  scale_color_manual(values=c("turquoise3")) +
  annotate("text", x = Inf, y = Inf, label = "N=2418", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item145: How much time in a normal day does DOG sped false digging (scratching at floor/carpet/etc)?", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 3, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 5), fill = "white"))


item146.hist1 <- ggplot(CCDitem146_fam, aes(x=item146)) +
  geom_histogram(alpha = 1,binwidth=1, color="thistle3", fill="thistle3") +
  scale_fill_manual(values=c("thistle3")) +
  scale_color_manual(values=c("thistle3")) +
  annotate("text", x = Inf, y = Inf, label = "N=2419", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item146: How much time in a normal day does DOG spend chasing HIS tail, or circling?", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())


table(CCDitem147_fam$item147, exclude = NULL)
item147.hist1 <- ggplot(CCDitem147_fam, aes(x=item147)) +
  geom_histogram(alpha = 1,binwidth=1, color="dodgerblue4", fill="dodgerblue4") +
  scale_fill_manual(values=c("dodgerblue4")) +
  scale_color_manual(values=c("dodgerblue4")) +
  annotate("text", x = Inf, y = Inf, label = "N=2414", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6, limits =c(0.5,6,5)) +
  ylim(0,2300) +
  labs(title = "item147: How much time in a normal day does DOG spend catching invisible flies?", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())



item148.hist1 <- ggplot(CCDitem148_fam, aes(x=item148)) +
  geom_histogram(alpha = 1,binwidth=1, color="darkseagreen3", fill="darkseagreen3") +
  scale_fill_manual(values=c("darkseagreen3")) +
  scale_color_manual(values=c("darkseagreen3")) +
  annotate("text", x = Inf, y = Inf, label = "N=2424", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item148: How much time in a normal day does DOG spend licking, chewing, or sucking on themself?", x = "", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 3, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 5), fill = "white"))

item149.hist1 <- ggplot(CCDitem149_fam, aes(x=item149)) +
  geom_histogram(alpha = 1,binwidth=1, color="darkslateblue", fill="darkslateblue") +
  scale_fill_manual(values=c("darkslateblue")) +
  scale_color_manual(values=c("darkslateblue")) +
  annotate("text", x = Inf, y = Inf, label = "N=2422", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item149: How much time in a normal day does DOG spend licking, chewing, or sucking on toys, people, or other pets?", 
       x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())

item151.hist1 <- ggplot(CCDitem151_fam, aes(x=item151)) +
  geom_histogram(alpha = 1,binwidth=1, color="coral1", fill="coral1") +
  scale_fill_manual(values=c("coral1")) +
  scale_color_manual(values=c("coral1")) +
  annotate("text", x = Inf, y = Inf, label = "N=2421", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item151: How much time in a normal day does DOG spend licking floors?", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 2), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 1, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 0), fill = "white"))+
  theme(axis.title.y=element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())



item152.hist1 <- ggplot(CCDitem152_fam, aes(x=item152)) +
  geom_histogram(alpha = 1,binwidth=1, color="dodgerblue2", fill="dodgerblue2") +
  scale_fill_manual(values=c("dodgerblue2")) +
  scale_color_manual(values=c("dodgerblue2")) +
  annotate("text", x = Inf, y = Inf, label = "N=2415", hjust = 1.3, vjust = 2.5, size = 4) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2300) +
  labs(title = "item152: How much time in a normal day does DOG spend chasing lights or reflections?", x = "Response", y = "Number of dogs")+
  theme(plot.title.position = "plot", 
        plot.title = element_textbox_simple(size = 10, padding = margin(5.5, 5.5, 5.5, 5.5), margin = margin(0, 0, 5.5, 42.5), 
                                            fill = "cornsilk"),
        axis.title.x = element_textbox_simple(width = NULL, padding = margin(4, 4, 4, 4), margin = margin(4, 0, 0, 0), 
                                              fill = "white"),
        axis.title.y = element_textbox_simple(hjust = 3, orientation = "left-rotated", minwidth = unit(1, "in"), maxwidth = unit(2, "in"),
                                              padding = margin(4, 4, 2, 4), margin = margin(0, 0, 2, 5), fill = "white"))


item145.hist1+item146.hist1+item147.hist1+item148.hist1+item149.hist1+item151.hist1+item152.hist1



