#Figure 3

#Geography

#PCoA Genus
#Genus and Corr PCoA Figure
#Extract eigen values (values of the variance reported) after running ordinate function with bray distance (used by phyloseq in generating PCoA plots)
b <- ordinate(Geography.G.tr, method="PCoA", dist="bray")
b.DF <- as.data.frame(b$vectors)
e.DF <- as.data.frame(b$values$Relative_eig)
e.DF["Percent"] <- (b$values$Relative_eig*100)

#Find out how many Axis needed to reach >= 1 or 100% variance
esum = 0
sum100 = 0
for (i in 1:nrow(e.DF))
{
  esum = esum + e.DF[i,1]
  if(esum >= 1)
  {
    sum100 = i
    break
  }
}

G <- Geography.G.tr
G.tr.DF <- as.data.frame(t(G@otu_table))

n = ncol(G.tr.DF)

G.tr.DF["Status"] <- NA

#Add the needed Axis Values

for(i in 1:sum100)
{
  G.tr.DF[colnames(b.DF[i])] <- NA
          for (j in 1:nrow(G.tr.DF))
          {
            G.tr.DF[j,colnames(b.DF[i])] <- b.DF[rownames(G.tr.DF[j,]),i]
          }
}

G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Kinshasa", "Masimanimba","Kahemba_Control_NonIntervention", "Kahemba_Control_Intervention"))

#G.tr.DF now has the rel abund data, with eigen values for each axis as additional columns so correlation can be done between genus relative abundance and PCoA Axis values to see which genus correlated the best with each Axis values
                                             
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Geography.G.tr@sam_data[rownames(G.tr.DF[i,]),]$Status
  }

a = n+2 #moves starting position to Axis.1 column
#for loop to correlate each species with each Axis.1 using the spearman method                                                                                          
for (i in a:ncol(G.tr.DF))
{
  Cor <- matrix(nrow = n,  ncol = 3)
  colnames(Cor) <- c(paste(colnames(G.tr.DF[i])), "spearman cor", "p-value")
  for (j in 1:n)
  {
    cor <- cor.test(G.tr.DF[,j], G.tr.DF[,i], method=c("spearman"))
    Cor[j,1] = colnames(G.tr.DF[j])
    Cor[j,2] = as.numeric(cor$estimate)
    Cor[j,3] = as.numeric(cor$p.value)
  }
  write.csv(Cor, file = paste("Geography_Genus_",colnames(G.tr.DF[i]),"_Correlation.csv", sep = ""))
  Cor <- data.frame(Cor, row.names = TRUE)
  Cor.f <- subset(Cor, rownames(Cor) %in% f_0.0001)                                       
  write.csv(Cor.f, file = paste("Geography_Genus_f_0.0001_",colnames(G.tr.DF[i]),"_Correlation.csv", sep = ""))
}

#Plot most correlated with Axis 1 and Axis 2                                             
#Correlation Plot
#Axis 1 Prevotella
#Axis 2 Lachnoclostridium

a1 <- ggplot(G.tr.DF, aes(x = Axis.1, y = Prevotella)) +
    geom_point(aes(color = factor(Status)), size = 1, stroke = 0, shape = 16) + theme_classic() + ylab("Prev.") + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 7))
a1 <- a1 + scale_color_manual(labels = SL, values = geography_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.y = element_text(size = 7, face = "italic"), axis.text.y = element_text(size = 7))  
a1 <- a1 + geom_smooth(method=lm, color = "black", size = 0.5) + theme(plot.margin=unit(c(0.15,0.15,0.25,0.15), "lines")) + scale_y_continuous(breaks = seq(-0.1, 0.4, by = 0.2))
#a1 <- a1 + stat_cor(method = "spearman", size = 5) 

a2 <- ggplot(G.tr.DF, aes(x = Lachnoclostridium, y = Axis.2)) +
    geom_point(aes(color = factor(Status)), size = 1, stroke = 0, shape = 16) + theme_classic() + xlab("Lach.") + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 7))
a2 <- a2 + scale_color_manual(labels = SL, values = geography_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.x = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  
a2 <- a2 + scale_y_continuous(position = "right") + scale_x_continuous(position = "top", breaks = seq(0.01, 0.02, by = 0.01))
a2 <- a2 + geom_smooth(method=lm, color = "black", size = 0.5) + theme(plot.margin=unit(c(0.15,0.125,0.15,0.15), "lines"))

p1 = plot_ordination(Geography.G.tr, ordinate(Geography.G.tr, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 2, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PGB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = geography_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA), legend.margin=margin(c(0,0,0,0))) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

l <- get_legend(PGB)
l <- as_ggplot(l)
l <- l + theme(plot.margin=unit(c(-1,0,0,-1), "lines"))

PGBt <- PGB + stat_ellipse(type = "t") + scale_x_continuous(position = "top") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
PGBt <- PGBt + theme(legend.position="none")
PGBt <- PGBt + annotate("text", x = -0.43, y = -0.42, label = expression(paste("p = 9.999x",10^-5)), size = 2.5)
PGBt <- ggarrange(PGBt,labels = c("A"),font.label = list(size = 7))

PGB <- PGB + scale_x_continuous(position = "top") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines")) + theme(legend.position="none")


G <- arrangeGrob(PGBt, a1,                               # bar plot spaning two columns
             a2, l,                               # box plot and scatter plot
             ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1,1,1,3), c(1,1,1,3), c(1,1,1,3), c(2, 2, 2, 4)))

tiff(filename = "Overall_Geography_Genus_PCoA_Corr.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
ggarrange(as_ggplot(G))
dev.off()


#PCoA KO
