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
#Creating phyloseq with only one group                                                 
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
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Class")

#OTU
Konzo_otu_c <- read.csv("./KinshasaControl_Konzo3_Bacteria_Class_ReadCounts.csv")
Konzo_class <- read.csv("./KinshasaControl_Konzo3_Bacteria_class.csv")
Konzo_Otu_C <-as.matrix(unname(Konzo_otu_c[1:nrow(Konzo_otu_c),5:(ncol(Konzo_otu_c))]))
rownames(Konzo_Otu_C)<-as.character( Konzo_otu_c[,1])
nam <-names(Konzo_otu_c)
colnames(Konzo_Otu_C)<-c(as.character(nam[5:length(nam)]))
OTU_C = otu_table(Konzo_Otu_C, taxa_are_rows = TRUE)

#TAX
Konzo_Class<-as.matrix(unname(Konzo_class[,2]))
rownames(Konzo_Class)<-as.character(unname(Konzo_Class[,1]))
colnames(Konzo_Class)<-"class"
TAX_C = tax_table(Konzo_Class)

#PhyloseqObject
KonzoData_C <-phyloseq(OTU_C, TAX_C, META)
#set all Na's to 0
KonzoData_C@otu_table[is.na(KonzoData_C@otu_table)] <- 0
KonzoData.C <- tax_glom(KonzoData_C, taxrank = "class")
KonzoData.C@sam_data$Status <- factor(KonzoData.C@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#Read Counts to Relative Abundance
KonzoData.C.tr <- transform_sample_counts(KonzoData.C, function(x) x / sum(x))

KonzoData.C.tr.status <- merge_samples(KonzoData.C.tr, KonzoData.C.tr@sam_data$Status)
KonzoData.C.tr.status <- transform_sample_counts(KonzoData.C.tr.status, function(x) x / 30)
                                                                                    
Kinshasa.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa", KonzoData.C)
Kinshasa.C.tr <- transform_sample_counts(Kinshasa.C, function(x) x / sum(x))
Masimanimba.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba", KonzoData.C)
Masimanimba.C.tr <- transform_sample_counts(Masimanimba.C, function(x) x / sum(x))                                        
ULPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.C)
ULPZ.C.tr <- transform_sample_counts(ULPZ.C, function(x) x / sum(x))
KLPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.C)
KLPZ.C.tr <- transform_sample_counts(KLPZ.C, function(x) x / sum(x))
UHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.C)
UHPZ.C.tr <- transform_sample_counts(UHPZ.C, function(x) x / sum(x))
KHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.C)
KHPZ.C.tr <- transform_sample_counts(KHPZ.C, function(x) x / sum(x))
                                     
#Filtering where mean in >= 0.01% (1E-4)

Kinshasa.C.tr.f <- filter_taxa(Kinshasa.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.C.tr.f <- filter_taxa(Masimanimba.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.C.tr.f <- filter_taxa(ULPZ.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.C.tr.f <- filter_taxa(KLPZ.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.C.tr.f <- filter_taxa(UHPZ.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.C.tr.f <- filter_taxa(KHPZ.C.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.C.tr.f@tax_table,Masimanimba.C.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.C.tr.f@tax_table, KLPZ.C.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.C.tr.f@tax_table,KHPZ.C.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

                           
write.csv(filterList, file = "Kinshasa_Konzo3_Class_f_0.0001.csv")
                           
x <- read.csv("Kinshasa_Konzo3_Class_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData.C.f <- prune_taxa(f_0.0001, KonzoData.C)
KonzoData.C.tr.f <- prune_taxa(f_0.0001, KonzoData.C.tr)

KonzoData.C.tr.status.f <- prune_taxa(f_0.0001, KonzoData.C.tr.status)

write.csv((KonzoData.C@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Class_ReadCounts.csv")
write.csv((KonzoData.C.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Class_RelAbund.csv")
write.csv(t(KonzoData.C.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Class_Avg_RelAbund.csv")
                                                                                 
                           
#Bacteria Order
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Order")

#OTU
Konzo_otu_o <- read.csv("./KinshasaControl_Konzo3_Bacteria_Order_ReadCounts.csv")
Konzo_order <- read.csv("./KinshasaControl_Konzo3_Bacteria_order.csv")
Konzo_Otu_O <-as.matrix(unname(Konzo_otu_o[1:nrow(Konzo_otu_o),5:(ncol(Konzo_otu_o))]))
rownames(Konzo_Otu_O)<-as.character( Konzo_otu_o[,1])
nam <-names(Konzo_otu_o)
colnames(Konzo_Otu_O)<-c(as.character(nam[5:length(nam)]))
OTU_O = otu_table(Konzo_Otu_O, taxa_are_rows = TRUE)

#TAX
Konzo_Order<-as.matrix(unname(Konzo_order[,2]))
rownames(Konzo_Order)<-as.character(unname(Konzo_Order[,1]))
colnames(Konzo_Order)<-"order"
TAX_O = tax_table(Konzo_Order)

#PhyloseqObject
KonzoData_O <-phyloseq(OTU_O, TAX_O, META)
#set all Na's to 0
KonzoData_O@otu_table[is.na(KonzoData_O@otu_table)] <- 0
KonzoData.O <- tax_glom(KonzoData_O, taxrank = "order")
KonzoData.O@sam_data$Status <- factor(KonzoData.O@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#Reads Counts to Relative Abundance
KonzoData.O.tr <- transform_sample_counts(KonzoData.O, function(x) x / sum(x))

KonzoData.O.tr.status <- merge_samples(KonzoData.O.tr, KonzoData.O.tr@sam_data$Status, fun = mean)
KonzoData.O.tr.status <- transform_sample_counts(KonzoData.O.tr.status, function(x) x / 30)
                           
#Filter
                                                 
Kinshasa.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa", KonzoData.O)
Kinshasa.O.tr <- transform_sample_counts(Kinshasa.O, function(x) x / sum(x))
Masimanimba.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba", KonzoData.O)
Masimanimba.O.tr <- transform_sample_counts(Masimanimba.O, function(x) x / sum(x))                                        
ULPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.O)
ULPZ.O.tr <- transform_sample_counts(ULPZ.O, function(x) x / sum(x))
KLPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.O)
KLPZ.O.tr <- transform_sample_counts(KLPZ.O, function(x) x / sum(x))
UHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.O)
UHPZ.O.tr <- transform_sample_counts(UHPZ.O, function(x) x / sum(x))
KHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.O)
KHPZ.O.tr <- transform_sample_counts(KHPZ.O, function(x) x / sum(x))
                                     
Kinshasa.O.tr.f <- filter_taxa(Kinshasa.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.O.tr.f <- filter_taxa(Masimanimba.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.O.tr.f <- filter_taxa(ULPZ.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.O.tr.f <- filter_taxa(KLPZ.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.O.tr.f <- filter_taxa(UHPZ.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.O.tr.f <- filter_taxa(KHPZ.O.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.O.tr.f@tax_table,Masimanimba.O.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.O.tr.f@tax_table, KLPZ.O.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.O.tr.f@tax_table,KHPZ.O.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(filterList, file = "Kinshasa_Konzo3_Order_f_0.0001.csv")

x <- read.csv("Kinshasa_Konzo3_Order_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
 
KonzoData.O.f <- prune_taxa(f_0.0001, KonzoData.O)
KonzoData.O.tr.f <- prune_taxa(f_0.0001, KonzoData.O.tr)
                                                 
KonzoData.O.tr.status.f <- prune_taxa(f_0.0001, KonzoData.O.tr.status)

write.csv((KonzoData.O@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Order_ReadCounts.csv")
write.csv((KonzoData.O.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Order_RelAbund.csv")
write.csv(t(KonzoData.O.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Order_Avg_RelAbund.csv")
                           
                           
#Bacteria Family
#FAMILY
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Family")

#OTU
Konzo_otu_f <- read.csv("./KinshasaControl_Konzo3_Bacteria_Family_ReadCounts.csv")
Konzo_family <- read.csv("./KinshasaControl_Konzo3_Bacteria_family.csv")
Konzo_Otu_F <-as.matrix(unname(Konzo_otu_f[1:nrow(Konzo_otu_f),5:(ncol(Konzo_otu_f))]))
rownames(Konzo_Otu_F)<-as.character( Konzo_otu_f[,1])
nam <-names(Konzo_otu_f)
colnames(Konzo_Otu_F)<-c(as.character(nam[5:length(nam)]))
OTU_F = otu_table(Konzo_Otu_F, taxa_are_rows = TRUE)

#TAX
Konzo_Family<-as.matrix(unname(Konzo_family[,2]))
rownames(Konzo_Family)<-as.character(unname(Konzo_Family[,1]))
colnames(Konzo_Family)<-"family"
TAX_F = tax_table(Konzo_Family)

#PhyloseqObject
KonzoData_F <-phyloseq(OTU_F, TAX_F, META)
#set all Na's to 0
KonzoData_F@otu_table[is.na(KonzoData_F@otu_table)] <- 0
KonzoData.F <- tax_glom(KonzoData_F, taxrank = "family")
KonzoData.F@sam_data$Status <- factor(KonzoData.F@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#Read Counts to Relative Abundance
KonzoData.F.tr <- transform_sample_counts(KonzoData.F, function(x) x / sum(x))

KonzoData.F.tr.status <- merge_samples(KonzoData.F.tr, KonzoData.F.tr@sam_data$Status)
KonzoData.F.tr.status <- transform_sample_counts(KonzoData.F.tr.status, function(x) x / 30)                                          
write.csv(t(KonzoData.F.tr.status@otu_table), file = "./KonzoDataFamily_AvgRelAbund_ByStatus.csv")
 
#Filter
                                                                                                                        
Kinshasa.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa", KonzoData.F)
Kinshasa.F.tr <- transform_sample_counts(Kinshasa.F, function(x) x / sum(x))
Masimanimba.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba", KonzoData.F)
Masimanimba.F.tr <- transform_sample_counts(Masimanimba.F, function(x) x / sum(x))                                        
ULPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.F)
ULPZ.F.tr <- transform_sample_counts(ULPZ.F, function(x) x / sum(x))
KLPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.F)
KLPZ.F.tr <- transform_sample_counts(KLPZ.F, function(x) x / sum(x))
UHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.F)
UHPZ.F.tr <- transform_sample_counts(UHPZ.F, function(x) x / sum(x))
KHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.F)
KHPZ.F.tr <- transform_sample_counts(KHPZ.F, function(x) x / sum(x))
                          
Kinshasa.F.tr.f <- filter_taxa(Kinshasa.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.F.tr.f <- filter_taxa(Masimanimba.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.F.tr.f <- filter_taxa(ULPZ.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.F.tr.f <- filter_taxa(KLPZ.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.F.tr.f <- filter_taxa(UHPZ.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.F.tr.f <- filter_taxa(KHPZ.F.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.F.tr.f@tax_table,Masimanimba.F.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.F.tr.f@tax_table, KLPZ.F.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.F.tr.f@tax_table,KHPZ.F.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(filterList, file = "Kinshasa_Konzo3_Family_f_0.0001.csv")
                                 
x <- read.csv("Kinshasa_Konzo3_Family_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.F.f <- prune_taxa(f_0.0001, KonzoData.F)
KonzoData.F.tr.f <- prune_taxa(f_0.0001, KonzoData.F.tr)

KonzoData.F.tr.status.f <- prune_taxa(f_0.0001, KonzoData.F.tr.status)

#write.csv((KonzoData.F@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Family_ReadCounts.csv")
#write.csv((KonzoData.F.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Family_RelAbund.csv")
#write.csv(t(KonzoData.F.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Family_Avg_RelAbund.csv")
                               
#Bacteria Genus
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")

#OTU
Konzo_otu_g <- read.csv("./KinshasaControl_Konzo3_Bacteria_Genus_ReadCounts.csv")
Konzo_genus <- read.csv("./KinshasaControl_Konzo3_Bacteria_genus.csv")
Konzo_Otu_G <-as.matrix(unname(Konzo_otu_g[1:nrow(Konzo_otu_g),5:(ncol(Konzo_otu_g))]))
rownames(Konzo_Otu_G)<-as.character( Konzo_otu_g[,1])
nam <-names(Konzo_otu_g)
colnames(Konzo_Otu_G)<-c(as.character(nam[5:length(nam)]))
OTU_G = otu_table(Konzo_Otu_G, taxa_are_rows = TRUE)

#TAX
Konzo_Genus<-as.matrix(unname(Konzo_genus[,2]))
rownames(Konzo_Genus)<-as.character(unname(Konzo_Genus[,1]))
colnames(Konzo_Genus)<-"genus"
TAX_G = tax_table(Konzo_Genus)

#PhyloseqObject
KonzoData_G <-phyloseq(OTU_G, TAX_G, META)

#set all Na's to 0
KonzoData_G@otu_table[is.na(KonzoData_G@otu_table)] <- 0

#Might not be necessary to taxglom since only genus in input file
KonzoData.G <- tax_glom(KonzoData_G, taxrank = "genus") 
KonzoData.G@sam_data$Status <- factor(KonzoData.G@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))
KonzoData.G@sam_data$Region <- factor(KonzoData.G@sam_data$Region, levels = c("Kinshasa", "Masimanimba", "Kahemba"))

#Read Count to Relative Abundance
KonzoData.G.tr <- transform_sample_counts(KonzoData.G, function(x) x / sum(x))
KonzoData.G.tr.log10 <- transform_sample_counts(KonzoData.G.tr, function(x) log10(x)) 


KonzoData.G.tr.status <- merge_samples(KonzoData.G.tr, KonzoData.G.tr@sam_data$Status, fun = mean)
KonzoData.G.tr.status <- transform_sample_counts(KonzoData.G.tr.status, function(x) x / 30)
                                                                       
#Filter                           
Kinshasa.G <- prune_samples(KonzoData.G@sam_data$Status == "Kinshasa", KonzoData.G)
Kinshasa.G.tr <- transform_sample_counts(Kinshasa.G, function(x) x / sum(x))
Masimanimba.G <- prune_samples(KonzoData.G@sam_data$Status == "Masimanimba", KonzoData.G)
Masimanimba.G.tr <- transform_sample_counts(Masimanimba.G, function(x) x / sum(x))                                        
ULPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.G)
ULPZ.G.tr <- transform_sample_counts(ULPZ.G, function(x) x / sum(x))
KLPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.G)
KLPZ.G.tr <- transform_sample_counts(KLPZ.G, function(x) x / sum(x))
UHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.G)
UHPZ.G.tr <- transform_sample_counts(UHPZ.G, function(x) x / sum(x))
KHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.G)
KHPZ.G.tr <- transform_sample_counts(KHPZ.G, function(x) x / sum(x))                           
                           
                          
Kinshasa.G.tr.f <- filter_taxa(Kinshasa.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.G.tr.f <- filter_taxa(Masimanimba.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.G.tr.f <- filter_taxa(ULPZ.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.G.tr.f <- filter_taxa(KLPZ.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.G.tr.f <- filter_taxa(UHPZ.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.G.tr.f <- filter_taxa(KHPZ.G.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.G.tr.f@tax_table,Masimanimba.G.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.G.tr.f@tax_table, KLPZ.G.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.G.tr.f@tax_table,KHPZ.G.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ


write.csv(filterList, file = "Kinshasa_Konzo3_Genus_f_0.0001.csv")

x <- read.csv("Kinshasa_Konzo3_Genus_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.G.f <- prune_taxa(f_0.0001, KonzoData.G)
KonzoData.G.tr.f <- prune_taxa(f_0.0001, KonzoData.G.tr)

KonzoData.G.tr.status.f <- prune_taxa(f_0.0001, KonzoData.G.tr.status)

write.csv((KonzoData.G@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Genus_ReadCounts.csv")
write.csv((KonzoData.G.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Genus_RelAbund.csv")
write.csv(t(KonzoData.G.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Genus_Avg_RelAbund.csv")
                           
                           
#Bacteria Species
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Species")

#OTU
Konzo_otu_s <- read.csv("./KinshasaControl_Konzo3_Bacteria_Species.ReadCounts.csv")
Konzo_species <- read.csv("./Kinshasa_Konzo3_species.csv")
Konzo_Otu_S <-as.matrix(unname(Konzo_otu_s[1:nrow(Konzo_otu_s),5:(ncol(Konzo_otu_s))]))

#Konzo_Otu_T_S <- t(Konzo_Otu_S)

rownames(Konzo_Otu_S)<-as.character( Konzo_otu_s[,1])
nam <-names(Konzo_otu_s)
colnames(Konzo_Otu_S)<-c(as.character(nam[5:length(nam)]))
OTU_S = otu_table(Konzo_Otu_S, taxa_are_rows = TRUE)

#TAX
Konzo_Species<-as.matrix(unname(Konzo_species[,2]))
rownames(Konzo_Species)<-as.character(unname(Konzo_Species[,1]))
colnames(Konzo_Species)<-"species"
TAX_S = tax_table(Konzo_Species)

#PhyloseqObject
KonzoData_S <-phyloseq(OTU_S, TAX_S, META)
#Set NAs to 0
KonzoData_S@otu_table[is.na(KonzoData_S@otu_table)] <- 0

#Probs can skip next step since the input file only contains species data so not require to grab/glom species taxa only
KonzoData.S <- tax_glom(KonzoData_S, taxrank = "species")

KonzoData.S@sam_data$Status <- factor(KonzoData.S@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Kahemba_Control_NonIntervention", "Kahemba_Konzo_NonIntervention", "Kahemba_Control_Intervention", "Kahemba_Konzo_Intervention"))

#Read Counts to Relative Abundance
KonzoData.S.tr <- transform_sample_counts(KonzoData.S, function(x) x / sum(x))
#Log base 10 Transform the rel abund values                                          
KonzoData.S.tr.log10 <- transform_sample_counts(KonzoData.S.tr, function(x) log10(x)) 
                                                
KonzoData.S.tr.status <- merge_samples(KonzoData.S.tr, KonzoData.S.tr@sam_data$Status)
KonzoData.S.tr.status <- transform_sample_counts(KonzoData.S.tr.status, function(x) x / 30)
                                         
#Filter                           
Kinshasa.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa", KonzoData.S)
Kinshasa.S.tr <- transform_sample_counts(Kinshasa.S, function(x) x / sum(x))
Masimanimba.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba", KonzoData.S)
Masimanimba.S.tr <- transform_sample_counts(Masimanimba.S, function(x) x / sum(x))                                        
ULPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kahemba_Control_NonIntervention", KonzoData.S)
ULPZ.S.tr <- transform_sample_counts(ULPZ.S, function(x) x / sum(x))
KLPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kahemba_Konzo_NonIntervention", KonzoData.S)
KLPZ.S.tr <- transform_sample_counts(KLPZ.S, function(x) x / sum(x))
UHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kahemba_Control_Intervention", KonzoData.S)
UHPZ.S.tr <- transform_sample_counts(UHPZ.S, function(x) x / sum(x))
KHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kahemba_Konzo_Intervention", KonzoData.S)
KHPZ.S.tr <- transform_sample_counts(KHPZ.S, function(x) x / sum(x))

                                         
Kinshasa.S.tr.f <- filter_taxa(Kinshasa.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.S.tr.f <- filter_taxa(Masimanimba.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.S.tr.f <- filter_taxa(ULPZ.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.S.tr.f <- filter_taxa(KLPZ.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.S.tr.f <- filter_taxa(UHPZ.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.S.tr.f <- filter_taxa(KHPZ.S.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.S.tr.f@tax_table,Masimanimba.S.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.S.tr.f@tax_table, KLPZ.S.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.S.tr.f@tax_table,KHPZ.S.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(filterList, file = "Kinshasa_Konzo3_Species_f_0.0001.csv")

x <- read.csv("Kinshasa_Konzo3_Species_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.S.f <- prune_taxa(f_0.0001, KonzoData.S)
KonzoData.S.tr.f <- prune_taxa(f_0.0001, KonzoData.S.tr)

KonzoData.S.tr.status.f <- prune_taxa(f_0.0001, KonzoData.S.tr.status)
                                                                                                       
write.csv((KonzoData.S.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Species_RelAbund.csv")
write.csv((KonzoData.S@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Species_ReadCounts.csv")
write.csv(t(KonzoData.S.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Species_Avg_RelAbund.csv")
                           
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





