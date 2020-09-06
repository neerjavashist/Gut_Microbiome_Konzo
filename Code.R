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

### Konzo Meta Data
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

#PhyloseqObject
KonzoData_P <-phyloseq(OTU_P, TAX_P, META)
#Set NAs to 0
KonzoData_P@otu_table[is.na(KonzoData_P@otu_table)] <- 0
KonzoData.P <- tax_glom(KonzoData_P, taxrank = "phylum")
KonzoData.P@sam_data$Status <- factor(KonzoData.P@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#Read Counts to Relative Abundance
KonzoData.P.tr <- transform_sample_counts(KonzoData.P, function(x) x / sum(x))
                                          
KonzoData.P.tr.status <- merge_samples(KonzoData.P.tr, KonzoData.P.tr@sam_data$Status)
KonzoData.P.tr.status <- transform_sample_counts(KonzoData.P.tr.status, function(x) x / 30) #average the sum of relabund in each group
  
#keep >= 0.01% in atleast one group
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

#write.csv(filterList, file = "Kinshasa_Konzo3_Phylum_f_0.0001.csv")
                                                                                                 
x <- read.csv("Kinshasa_Konzo3_Phylum_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData.P.f <- prune_taxa(f_0.0001, KonzoData.P)
KonzoData.P.tr.f <- prune_taxa(f_0.0001, KonzoData.P.tr)                                                 
KonzoData.P.tr.status.f <- prune_taxa(f_0.0001, KonzoData.P.tr.status)
                                                 

#Bacteria Class
#Bacteria Order
#Bacteria Family
#Bacteria Genus
#Bacteria Species

### Estimate Richness

### Beta Diversity using Bray-Curtis for Bacteria Genus

### Mann Whitney-Wilcox Test (with BH correction)
#Bacteria Phylum
#Bacteria Class
#Bacteria Order
#Bacteria Family
#Bacteria Genus
#Bacteria Species





