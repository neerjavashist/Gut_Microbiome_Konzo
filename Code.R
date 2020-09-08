### All packages that may be required
library("patchwork")
library("ggpubr")
library("devtools")
library("picante")
library("TSA")
library("nortest")
library("multcomp")
library("car")
library("microbiome")
library("mvabund")
library("MASS")
library("geepack")
library("doBy")
library("lattice")
library("MuMIn")
library("phangorn")
library("DESeq2")
library("FSA")
library("phyloseq")
library("microbiome")
library("tidyverse")
library("ggplot2")
library("rentrez")
library("XML")
library("ggrepel")
library("plyr")
library("vegan")
library("DESeq2")
library("reshape2")
library("colorspace")
library("zoo")
library("magrittr")
library("tidyr")
library("cvcqv")
library("raster")

### The general format for each taxanomic rank is the same, so the code is repetitive. Comments are provided in the Bacteria Phylum section to show what is happening, and the explainations are applicable to the other taxa ranks as well.

### Naming Scheme Shorthand

SL <- c(Kinshasa = "Kinshasa", Masimanimba = "Masi-manimba", Kahemba_Control_NonIntervention = "Unaffected LPZ", Kahemba_Konzo_NonIntervention = "Konzo LPZ", Kahemba_Control_Intervention = "Unaffected HPZ", Kahemba_Konzo_Intervention = "Konzo HPZ")

SSL <- c(Kinshasa = "Kin", Masimanimba = "Mas", Kahemba_Control_NonIntervention = "Unaffected LPZ", Kahemba_Konzo_NonIntervention = "Konzo LPZ", Kahemba_Control_Intervention = "Unaffected HPZ", Kahemba_Konzo_Intervention = "Konzo HPZ")

SSSL <- c(Kinshasa = "Kin", Masimanimba = "Mas", Kahemba_Control_NonIntervention = "ULPZ", Kahemba_Konzo_NonIntervention = "KLPZ", Kahemba_Control_Intervention = "UHPZ", Kahemba_Konzo_Intervention = "KHPZ")

### Color Scheme

#Kinshasa: "royalblue1"
#Masi-manimba: "springgreen3"
#ULPZ: "turquoise3"
#KLPZ: "tomato"
#UHPZ: "slateblue1"
#KHPZ: "gold"

kinmas_color <- c("royalblue1",   "springgreen3")
kinulpz_color <- c("royalblue1",    "turquoise3")
masulpz_color <- c("springgreen3",     "turquoise3")
kinuhpz_color <- c("royalblue1",  "slateblue1")
masuhpz_color <- c("springgreen3","slateblue1")

nonintervention_color <- c("turquoise3",        "tomato" )
intervention_color <- c( "slateblue1",        "gold")
control_color <- c( "turquoise3",        "slateblue1")
disease_color <- c("tomato", "gold")
geography_color <- c("royalblue1",   "springgreen3", "turquoise3", "slateblue1")
kahemba_color <- c("turquoise3", "tomato", "slateblue1", "gold")
konzo_color <- c("royalblue1",   "springgreen3", "turquoise3", "tomato", "slateblue1", "gold")

### Konzo Meta Data
#Konzo_meta contains any additional information needed. The relevant data for the project are in the Supplemental File 1 Sample_Metadata tab. Columns Sample, Name, Run, ID, Region, Age, and Sex are as here.
#Konzo_meta "Status" is the Sample_Metadata "Group", and values Kahemba_Control_NonIntervention is changed to Kahemba_Unaffected_LPZ, Kahemba_Konzo_NonIntervention is changed to Kahemba_Konzo_LPZ, Kahemba_Control_Intervention is changed to Kahemba_Unaffected_HPZ, and Kahemba_Konzo_Intervention is changed to Kahemba_Konzo_HPZ

setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken")

#META
Konzo_meta <- read.csv("./KinshasaControl_Konzo3_Meta.csv")
names(Konzo_meta)<-c("Sample","Name","Run","ID","Region","Status","Disease","Sample_ID","Collection_date","DNA_Concentration","Isolation_date","Elution","Age","Sex","Disease_Old","Intervention")
rownames(Konzo_meta)<-as.character(Konzo_meta[,1])
META<-sample_data(Konzo_meta)

### Read Count to Relative Abundance
#Bacteria Phylum
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Phylum")

#OTU
Konzo_otu_p <- read.csv("./KinshasaControl_Konzo3_Bacteria_Phylum.ReadCounts.csv")
#Konzo_phylum contains the names of all the taxa in the ReadCounts file in one column with an empty first entry. The taxa names are duplicated into the second column and the second column has a column name although this will be removed in the code. 
Konzo_phylum <- read.csv("./Kinshasa_Konzo3_phylum.csv")
Konzo_Otu_P <-as.matrix(unname(Konzo_otu_p[1:nrow(Konzo_otu_p),5:(ncol(Konzo_otu_p))]))
rownames(Konzo_Otu_P)<-as.character( Konzo_otu_p[,1])
nam <-names(Konzo_otu_p)
colnames(Konzo_Otu_P)<-c(as.character(nam[5:length(nam)]))
OTU_P = otu_table(Konzo_Otu_P, taxa_are_rows = TRUE)

#TAX
Konzo_Phylum<-as.matrix(unname(Konzo_phylum[,2]))
rownames(Konzo_Phylum)<-as.character(unname(Konzo_Phylum[,1]))
colnames(Konzo_Phylum)<-"phylum"
TAX_P = tax_table(Konzo_Phylum)

#Create the PhyloseqObject
KonzoData_P <-phyloseq(OTU_P, TAX_P, META)
#Set NAs to 0
KonzoData_P@otu_table[is.na(KonzoData_P@otu_table)] <- 0
KonzoData.P <- tax_glom(KonzoData_P, taxrank = "phylum")
KonzoData.P@sam_data$Status <- factor(KonzoData.P@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))
#Writing the otu_table. Supplemental File 1, Phylum Tab
write.csv(KonzoData.P@otu_table), file = "./KonzoDataPhylum_ReadCounts.csv")  


#Read Counts to Relative Abundance
KonzoData.P.tr <- transform_sample_counts(KonzoData.P, function(x) x / sum(x))
#Writing the otu_table. Supplemental File 2, Phylum Tab                                     
write.csv(KonzoData.P.tr@otu_table), file = "./KonzoDataPhylum_AvgRelAbund.csv")  
#Merge samples by group/status                                         
KonzoData.P.tr.status <- merge_samples(KonzoData.P.tr, KonzoData.P.tr@sam_data$Status) #merge_smaples by default sums the values for otu
KonzoData.P.tr.status <- transform_sample_counts(KonzoData.P.tr.status, function(x) x / 30) #average the sum of relabund in each group
                                                 
#Writing the otu_table. Supplemental File 2, Phylum Tab                                                                                      
write.csv(t(KonzoData.P.tr.status@otu_table), file = "./KonzoDataPhylum_AvgRelAbund_ByStatus.csv")
  
#keep Rel abund >= 0.01% in atleast one group
Kinshasa.P <- prune_samples(KonzoData.P@sam_data$Status == "Kinshasa", KonzoData.P)
Kinshasa.P.tr <- transform_sample_counts(Kinshasa.P, function(x) x / sum(x))
Masimanimba.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba", KonzoData.P)
Masimanimba.P.tr <- transform_sample_counts(Masimanimba.P, function(x) x / sum(x))                                        
ULPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.P)
ULPZ.P.tr <- transform_sample_counts(ULPZ.P, function(x) x / sum(x))
KLPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.P)
KLPZ.P.tr <- transform_sample_counts(KLPZ.P, function(x) x / sum(x))
UHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.P)
UHPZ.P.tr <- transform_sample_counts(UHPZ.P, function(x) x / sum(x))
KHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.P)
KHPZ.P.tr <- transform_sample_counts(KHPZ.P, function(x) x / sum(x))
                                     
Kinshasa.P.tr.f <- filter_taxa(Kinshasa.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.P.tr.f <- filter_taxa(Masimanimba.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.P.tr.f <- filter_taxa(ULPZ.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.P.tr.f <- filter_taxa(KLPZ.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.P.tr.f <- filter_taxa(UHPZ.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.P.tr.f <- filter_taxa(KHPZ.P.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.P.tr.f@tax_table,Masimanimba.P.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.P.tr.f@tax_table, KLPZ.P.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.P.tr.f@tax_table,KHPZ.P.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ
#Save the filter list for future filtering
write.csv(filterList, file = "Kinshasa_Konzo3_Phylum_f_0.0001.csv")
                                                                                                 
x <- read.csv("Kinshasa_Konzo3_Phylum_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData.P.f <- prune_taxa(f_0.0001, KonzoData.P) #filtered readcount phyloseq object
KonzoData.P.tr.f <- prune_taxa(f_0.0001, KonzoData.P.tr) #filtered rel abund phyloseq object                                            
KonzoData.P.tr.status.f <- prune_taxa(f_0.0001, KonzoData.P.tr.status) #filtered rel abund megerd by groups/status phyoseq object
                                                 
#Bacteria Class
#Bacteria Order
#Bacteria Family
#Bacteria Genus
#Bacteria Species

### Estimate Richness
#Read Count from KonzoData.S (Bacteria Species data)                           
otuD.S <- as.data.frame(t(otu_table(KonzoData.S)))
diversity.S <- estimate_richness(KonzoData.S)
diversity.S <- cbind(sample_data(KonzoData.S),diversity.S) #Check if correct sample data was cbind. Can be tricky so always confirm
diversity.S$Status <- as.factor(diversity.S$Status)
diversity.S$Status <- factor(diversity.S$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#STATISTICS for Estimate Richness
#One-way ANOVA to see if there is a statitically significant difference in the measure of alpha diversity
observed.aov <- aov(Observed ~ Status, data = diversity.S)
chao1.aov <- aov(Chao1 ~ Status, data = diversity.S)
shannon.aov <- aov(Shannon ~ Status, data = diversity.S)
ACE.aov <- aov(ACE ~ Status, data = diversity.S)
simpson.aov <- aov(Simpson ~ Status, data = diversity.S)
fisher.aov <- aov(Fisher ~ Status, data = diversity.S)
   
###Figure 2:                           
observed <- ggplot(diversity.S, aes(factor(Status), Observed)) + geom_boxplot(aes(fill = factor(Status)),fatten = 1, outlier.shape = NA) + labs(x = element_blank(), y = "OTU") + theme(axis.text.x = element_blank()) + theme_classic()
observed <- observed + geom_jitter(position=position_jitter(0.2), size = 0.3)
observed2 <- observed + stat_summary(fun=mean, geom="point", shape=23, size=1.5, color = "black", fill="white")
observed3 <- observed2 + theme(legend.position="bottom", legend.margin=margin(-10,0,0,0)) + theme(legend.direction = "horizontal") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(axis.ticks.x = element_blank(), axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_blank())

observed4 <- observed3 + guides(fill=guide_legend(ncol=6)) 
observed4 <- observed4 + scale_fill_manual(values = konzo_color, labels = SSSL)

shan <- ggplot(diversity.S, aes(factor(Status), Shannon))+ geom_boxplot(aes(fill = factor(Status)),fatten = 1, outlier.shape = NA) + labs(x = element_blank(), y = "Shannon Diversity Index") + theme(axis.text.x = element_blank()) + theme_classic()
shan <- shan + geom_jitter(position=position_jitter(0.2), size = 0.3)
shan2 <- shan + stat_summary(fun=mean, geom="point", shape=23, size=1.5, color = "black", fill="white")
shan3 <- shan2 + theme(legend.position="bottom", legend.margin=margin(-10,0,0,0)) + theme(legend.direction = "horizontal") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(axis.ticks.x = element_blank(), axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_blank())

shan4 <- shan3 + guides(fill=guide_legend(ncol=6)) 
shan4 <- shan4 + scale_fill_manual(values = konzo_color, labels = SSSL)


#Phylum
top_P <- read.csv("Kinshasa_Konzo3_Phylum_Top4.csv", row.names = 1, colClasses = "character")
top_P <- unlist(top_P)

KonzoData.P.tr.status.top = prune_taxa(top_P, KonzoData.P.tr.status)

p <- plot_bar(KonzoData.P.tr.status.top, "Sample", "Abundance", fill = 'phylum')
p$data$Sample <- factor(p$data$Sample, levels = c("Kahemba_Konzo_Intervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_NonIntervention", "Masimanimba", "Kinshasa"))

p <- p + labs(x = element_blank(), y = "Relative Abundance") +  scale_fill_discrete(name = "Phylum")
top_phylum_plot <- p + theme(legend.position="bottom") + theme(legend.key.size = unit(.2, "cm"))

top_phylum_plot <- top_phylum_plot + 
  theme(legend.position="right", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-8,0,-8,-8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_x_discrete(labels= SSSL)+
  theme(plot.title = element_blank(), legend.key.size = unit(.2, "cm"), legend.text = element_text(size = 7),legend.title = element_text(size = 7)) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme(axis.text.y = element_text(angle = 0, size = 7), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 0,vjust=1, hjust=0.5, size = 7), axis.title.y = element_text(size = 7), legend.title = element_text(size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_text(size=7))
top_phylum_plot <- top_phylum_plot + geom_bar(stat = "identity") + coord_flip() + theme(axis.title.x = element_text(size = 7))
top_phylum_plot

#Family
top_F <- read.csv("Kinshasa_Konzo3_Family_Top5.csv", row.names = 1, colClasses = "character")
top_F <- unlist(top_F)

KonzoData.F.tr.status.top = prune_taxa(top_F, KonzoData.F.tr.status)

p <- plot_bar(KonzoData.F.tr.status.top, "Sample", "Abundance", fill = 'family')
p$data$Sample <- factor(p$data$Sample, levels = c("Kahemba_Konzo_Intervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_NonIntervention", "Masimanimba", "Kinshasa"))

p <- p + labs(x = element_blank(), y = "Relative Abundance") +  scale_fill_discrete(name = "Family")
top_family_plot <- p + theme(legend.position="bottom") + theme(legend.key.size = unit(.2, "cm"), legend.text = element_text(size = 7, face="italic"))

top_family_plot <- top_family_plot + 
  theme(legend.position="right", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-8,0,-8,-8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.7)) +
  scale_x_discrete(labels= SSSL)+
  theme(plot.title = element_blank(), legend.key.size = unit(.2, "cm"), legend.text = element_text(size = 7), legend.title = element_text(size = 7)) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme(axis.text.y = element_text(angle = 0, size = 7), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 0,vjust=1, hjust=0.5, size = 7), axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), legend.title = element_text(size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_text(size=7))
top_family_plot <- top_family_plot + geom_bar(stat = "identity") + coord_flip()
top_family_plot

#Genus
top_G <- read.csv("Kinshasa_Konzo3_Genus_Top7.csv", row.names = 1, colClasses = "character")
top_G <- unlist(top_G)

################ Needs to be added
                           
#When adding heat map in gimp to full figure, make sure Image > Print Size has correct inches and ppi (set to >=300))

o <- as.data.frame(otu_table(KonzoData.S.tr.status.f))                                                 
tiff(filename = "KinshasaKonzo3_Bacteria_Species_Heatmap.tiff", width = 2.5, height = 3.25, units = "in", res = 600)
heatmap.2(as.matrix(t(o)), scale = "row", trace = "none", keysize = 0.25, labRow = "Species", labCol = SSSL, margins = c(1, 1), Rowv = FALSE, dendrogram = "column", key.title = NA, srtCol = 0, srtRow = 90 , cexCol = 0.75, cexRow = 0.75, offsetRow = 0, offsetCol = 0, lhei = c(0.5,2,2,1.25), lwid = c(0.1,1,1), key.par = list(cex=0.5), lmat = rbind(c(0,3,3),c(2,1,1),c(2,1,1),c(0,0,4)), adjCol = c(0.5,0.5), adjRow = c(4.5,0.25))
dev.off()                          
                           
### Beta Diversity using Bray-Curtis for Bacteria Genus
                           
                           

### Mann Whitney-Wilcox Test (with BH correction)
                           
#Bacteria Phylum
#Bacteria Class
#Bacteria Order
#Bacteria Family
#Bacteria Genus
#Bacteria Species





