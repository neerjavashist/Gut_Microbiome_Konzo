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
library(gridExtra)
library(grid)
library(lattice)
library(gplots)


### for pavian


### The general format for each taxanomic rank is the same, so the code is repetitive. Comments are provided in the Bacteria Phylum section to show what is happening, and the explainations are applicable to the other taxa ranks as well.

### For figures, since the standard R ggplot cannot generate corrected p-values, the p-values post BH correction must be manually added to the plots.

### Due to the way the input data is set up, the ordering of the actual samples within excel sheets may differ, so always check which columns correspond to which samples to be sure the correst samples are being compared (ex: for aldex)


### Naming Scheme Shorthand

SL <- c(Kinshasa = "Kinshasa", Masimanimba = "Masi-Manimba", Unaffected_Low_Prevalence_Zone = "Unaffected LPZ", Konzo_Low_Prevalence_Zone = "Konzo LPZ", Unaffected_High_Prevalence_Zone = "Unaffected HPZ", Konzo_High_Prevalence_Zone = "Konzo HPZ")

SSL <- c(Kinshasa = "Kin", Masimanimba = "Mas", Unaffected_Low_Prevalence_Zone = "Unaffected LPZ", Konzo_Low_Prevalence_Zone = "Konzo LPZ", Unaffected_High_Prevalence_Zone = "Unaffected HPZ", Konzo_High_Prevalence_Zone = "Konzo HPZ")

SSSL <- c(Kinshasa = "Kin", Masimanimba = "Mas", Unaffected_Low_Prevalence_Zone = "ULPZ", Konzo_Low_Prevalence_Zone = "KLPZ", Unaffected_High_Prevalence_Zone = "UHPZ", Konzo_High_Prevalence_Zone = "KHPZ")

### Color Scheme

#Kinshasa: "royalblue1"
#Masi-manimba: "springgreen3"
#ULPZ: "gold"
#KLPZ: "tomato"
#UHPZ: "mediumorchid1"
#KHPZ: "deepskyblue1"


#COLOR SCHEME
kinmas_color <- c("royalblue1",   "springgreen3")
kinulpz_color <- c("royalblue1",    "gold")
masulpz_color <- c("springgreen3",     "gold")
kinuhpz_color <- c("royalblue1",  "mediumorchid1")
masuhpz_color <- c("springgreen3","mediumorchid1")
kinklpz_color <- c("royalblue1",    "tomato" )
masklpz_color <- c("springgreen3",    "tomato" )
kinkhpz_color <- c("royalblue1",    "deepskyblue1")
maskhpz_color <- c("springgreen3",     "deepskyblue1")


lpz_color <- c("gold",        "tomato" )
hpz_color <- c( "mediumorchid1",        "deepskyblue1")
control_color <- c( "gold",        "mediumorchid1")
disease_color <- c("tomato", "deepskyblue1")
geography_color <- c("royalblue1",   "springgreen3", "gold", "mediumorchid1")
kahemba_color <- c("gold", "tomato", "mediumorchid1", "deepskyblue1")
konzo_color <- c("royalblue1",   "springgreen3", "gold", "tomato", "mediumorchid1", "deepskyblue1")



### Konzo Meta Data
#Konzo_meta contains any additional information needed. The relevant data for the project are in the Supplemental File 1 Sample_Metadata tab. Columns Sample, Name, Run, ID, Region, Age, and Sex are as here.
#Konzo_meta "Status" is the Sample_Metadata "Group"
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken")

#META
Konzo_meta <- read.csv("./KinshasaControl_Konzo3_Meta_Mod.csv")
names(Konzo_meta)<-c("Sample","Name","Run","ID","Region","Status","Disease","Sample_ID","Collection_date","DNA_Concentration","Isolation_date","Elution","Age","Sex","Disease_Old","Geography")
rownames(Konzo_meta)<-as.character(Konzo_meta[,1])
META<-sample_data(Konzo_meta)

### Read Count to Relative Abundance
#Bacteria Phylum
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Phylum")

#OTU
Konzo_otu_p <- read.csv("./KinshasaControl_Konzo3_Bacteria_Phylum.ReadCounts.csv") #refer to this file in the repository for reference (has additional columns that are removed in the code because they are unnecessary)
#Konzo_phylum contains the names of all the taxa in the ReadCounts file in one column with an empty first entry. The taxa names are duplicated into the second column and the second column has a column name although this will be removed in the code. 
Konzo_phylum <- read.csv("./Kinshasa_Konzo3_phylum.csv") #refer to this file in the repository for reference
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
KonzoData.P@sam_data$Status <- factor(KonzoData.P@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
#Writing the otu_table. Supplemental File 1, Phylum Tab
#write.csv(KonzoData.P@otu_table), file = "./KonzoDataPhylum_ReadCounts.csv")  


#Read Counts to Relative Abundance
KonzoData.P.tr <- transform_sample_counts(KonzoData.P, function(x) x / sum(x))
#write.csv(KonzoData.P.tr@otu_table), file = "./KonzoDataPhylum_AvgRelAbund.csv")  
#Merge samples by group/status                                         
KonzoData.P.tr.status <- merge_samples(KonzoData.P.tr, KonzoData.P.tr@sam_data$Status) #merge_smaples by default sums the values for otu
KonzoData.P.tr.status <- transform_sample_counts(KonzoData.P.tr.status, function(x) x / 30) #average the sum of relabund in each group
                                                 
#write.csv(t(KonzoData.P.tr.status@otu_table), file = "./KonzoDataPhylum_AvgRelAbund_ByStatus.csv")
                                                 
#keep Rel abund >= 0.01% in atleast one group
#Creating phyloseq with only one group                                                 
Kinshasa.P <- prune_samples(KonzoData.P@sam_data$Status == "Kinshasa", KonzoData.P)
Kinshasa.P.tr <- transform_sample_counts(Kinshasa.P, function(x) x / sum(x))
Masimanimba.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba", KonzoData.P)
Masimanimba.P.tr <- transform_sample_counts(Masimanimba.P, function(x) x / sum(x))                                        
ULPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.P)
ULPZ.P.tr <- transform_sample_counts(ULPZ.P, function(x) x / sum(x))
KLPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.P)
KLPZ.P.tr <- transform_sample_counts(KLPZ.P, function(x) x / sum(x))
UHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.P)
UHPZ.P.tr <- transform_sample_counts(UHPZ.P, function(x) x / sum(x))
KHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
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
#write.csv(filterList, file = "Kinshasa_Konzo3_Phylum_f_0.0001.csv")
                                                                                                 
x <- read.csv("Kinshasa_Konzo3_Phylum_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData.P.f <- prune_taxa(f_0.0001, KonzoData.P) #filtered readcount phyloseq object
KonzoData.P.tr.f <- prune_taxa(f_0.0001, KonzoData.P.tr) #filtered rel abund phyloseq object                                            
KonzoData.P.tr.status.f <- prune_taxa(f_0.0001, KonzoData.P.tr.status) #filtered rel abund megerd by groups/status phyoseq object

#Mean and Standard Deviation
KonzoData.P.tr.df <- as.data.frame(t(KonzoData.P.tr@otu_table))
KonzoData.P.tr.df <- cbind(KonzoData.P.tr.df, KonzoData.P.tr@sam_data$Status)

colnames(KonzoData.P.tr.df)[colnames(KonzoData.P.tr.df)=="KonzoData.P.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.P.tr.df))
  {KonzoData.P.tr.df[i,]$Status <- KonzoData.P.tr@sam_data[rownames(KonzoData.P.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.P.tr.avg <- KonzoData.P.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.P.tr.avg.x <- t(KonzoData.P.tr.avg)
colnames(KonzoData.P.tr.avg.x) <- KonzoData.P.tr.avg.x[1,]   
KonzoData.P.tr.avg.x <- KonzoData.P.tr.avg.x[-1,]
colnames(KonzoData.P.tr.avg.x) <- paste("Avg", colnames(KonzoData.P.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.P.tr.sd <- KonzoData.P.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.P.tr.sd.x <- t(KonzoData.P.tr.sd)
colnames(KonzoData.P.tr.sd.x) <- KonzoData.P.tr.sd.x[1,]   
KonzoData.P.tr.sd.x <- KonzoData.P.tr.sd.x[-1,]
colnames(KonzoData.P.tr.sd.x) <- paste("SD", colnames(KonzoData.P.tr.sd.x), sep = "_")                                                 

KonzoData.P.tr.avg.sd <-merge(KonzoData.P.tr.avg.x,KonzoData.P.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.P.tr.avg.sd) <- KonzoData.P.tr.avg.sd[,1]   
KonzoData.P.tr.avg.sd <- KonzoData.P.tr.avg.sd[,-1]  
                           
KonzoData.P.tr.avg.sd <- KonzoData.P.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.P.tr.avg.sd, file = "./KonzoDataPhylum_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.P.tr.avg.sd.f <- subset(KonzoData.P.tr.avg.sd, rownames(KonzoData.P.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.P.tr.avg.sd.f, file = "./KonzoDataPhylum_AvgRelAbund_SD_ByGroup_filtered.csv")                            

#Data for Relative Abundnace for Phylum (Supplemental File 2)                           
Phylum.tr <- merge(KonzoData.P.tr.avg.sd,as.data.frame(KonzoData.P.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Phylum.tr, file = "./KonzoDataPhylum_RelAbund_Supp.csv")                           
#Data for Relative Abundnace for Phylum Filtered (Supplemental File 2) 
Phylum.tr.f <- merge(KonzoData.P.tr.avg.sd.f,as.data.frame(KonzoData.P.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Phylum.tr.f, file = "./KonzoDataPhylum_RelAbund_filtered_Supp.csv")                           
                           
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
KonzoData.C@sam_data$Status <- factor(KonzoData.C@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
#Read Counts to Relative Abundance
KonzoData.C.tr <- transform_sample_counts(KonzoData.C, function(x) x / sum(x))

KonzoData.C.tr.status <- merge_samples(KonzoData.C.tr, KonzoData.C.tr@sam_data$Status)
KonzoData.C.tr.status <- transform_sample_counts(KonzoData.C.tr.status, function(x) x / 30)
                                                                                                 
                                                 
Kinshasa.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa", KonzoData.C)
Kinshasa.C.tr <- transform_sample_counts(Kinshasa.C, function(x) x / sum(x))
Masimanimba.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba", KonzoData.C)
Masimanimba.C.tr <- transform_sample_counts(Masimanimba.C, function(x) x / sum(x))                                        
ULPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.C)
ULPZ.C.tr <- transform_sample_counts(ULPZ.C, function(x) x / sum(x))
KLPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.C)
KLPZ.C.tr <- transform_sample_counts(KLPZ.C, function(x) x / sum(x))
UHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.C)
UHPZ.C.tr <- transform_sample_counts(UHPZ.C, function(x) x / sum(x))
KHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)
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

                           
#write.csv(filterList, file = "Kinshasa_Konzo3_Class_f_0.0001.csv")
                           
x <- read.csv("Kinshasa_Konzo3_Class_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData.C.f <- prune_taxa(f_0.0001, KonzoData.C)
KonzoData.C.tr.f <- prune_taxa(f_0.0001, KonzoData.C.tr)

KonzoData.C.tr.status.f <- prune_taxa(f_0.0001, KonzoData.C.tr.status)

#write.csv((KonzoData.C@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Class_ReadCounts.csv")
#write.csv((KonzoData.C.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Class_RelAbund.csv")
#write.csv(t(KonzoData.C.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Class_Avg_RelAbund.csv")

#Mean and Standard Deviation
KonzoData.C.tr.df <- as.data.frame(t(KonzoData.C.tr@otu_table))
KonzoData.C.tr.df <- cbind(KonzoData.C.tr.df, KonzoData.C.tr@sam_data$Status)

colnames(KonzoData.C.tr.df)[colnames(KonzoData.C.tr.df)=="KonzoData.C.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.C.tr.df))
  {KonzoData.C.tr.df[i,]$Status <- KonzoData.C.tr@sam_data[rownames(KonzoData.C.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.C.tr.avg <- KonzoData.C.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.C.tr.avg.x <- t(KonzoData.C.tr.avg)
colnames(KonzoData.C.tr.avg.x) <- KonzoData.C.tr.avg.x[1,]   
KonzoData.C.tr.avg.x <- KonzoData.C.tr.avg.x[-1,]
colnames(KonzoData.C.tr.avg.x) <- paste("Avg", colnames(KonzoData.C.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.C.tr.sd <- KonzoData.C.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.C.tr.sd.x <- t(KonzoData.C.tr.sd)
colnames(KonzoData.C.tr.sd.x) <- KonzoData.C.tr.sd.x[1,]   
KonzoData.C.tr.sd.x <- KonzoData.C.tr.sd.x[-1,]
colnames(KonzoData.C.tr.sd.x) <- paste("SD", colnames(KonzoData.C.tr.sd.x), sep = "_")                                                 

KonzoData.C.tr.avg.sd <-merge(KonzoData.C.tr.avg.x,KonzoData.C.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.C.tr.avg.sd) <- KonzoData.C.tr.avg.sd[,1]   
KonzoData.C.tr.avg.sd <- KonzoData.C.tr.avg.sd[,-1]  
                           
KonzoData.C.tr.avg.sd <- KonzoData.C.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.C.tr.avg.sd, file = "./KonzoDataClass_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.C.tr.avg.sd.f <- subset(KonzoData.C.tr.avg.sd, rownames(KonzoData.C.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.C.tr.avg.sd.f, file = "./KonzoDataClass_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
Class.tr <- merge(KonzoData.C.tr.avg.sd,as.data.frame(KonzoData.C.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Class.tr, file = "./KonzoDataClass_RelAbund_Supp.csv")                           

Class.tr.f <- merge(KonzoData.C.tr.avg.sd.f,as.data.frame(KonzoData.C.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Class.tr.f, file = "./KonzoDataClass_RelAbund_filtered_Supp.csv")                           
                           
                           
                           
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
KonzoData.O@sam_data$Status <- factor(KonzoData.O@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
#Reads Counts to Relative Abundance
KonzoData.O.tr <- transform_sample_counts(KonzoData.O, function(x) x / sum(x))

KonzoData.O.tr.status <- merge_samples(KonzoData.O.tr, KonzoData.O.tr@sam_data$Status, fun = mean)
KonzoData.O.tr.status <- transform_sample_counts(KonzoData.O.tr.status, function(x) x / 30)
                                                                                                                                                
#Filter
                                                 
Kinshasa.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa", KonzoData.O)
Kinshasa.O.tr <- transform_sample_counts(Kinshasa.O, function(x) x / sum(x))
Masimanimba.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba", KonzoData.O)
Masimanimba.O.tr <- transform_sample_counts(Masimanimba.O, function(x) x / sum(x))                                        
ULPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.O)
ULPZ.O.tr <- transform_sample_counts(ULPZ.O, function(x) x / sum(x))
KLPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.O)
KLPZ.O.tr <- transform_sample_counts(KLPZ.O, function(x) x / sum(x))
UHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.O)
UHPZ.O.tr <- transform_sample_counts(UHPZ.O, function(x) x / sum(x))
KHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)
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

#write.csv(filterList, file = "Kinshasa_Konzo3_Order_f_0.0001.csv")

x <- read.csv("Kinshasa_Konzo3_Order_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
 
KonzoData.O.f <- prune_taxa(f_0.0001, KonzoData.O)
KonzoData.O.tr.f <- prune_taxa(f_0.0001, KonzoData.O.tr)
                                                 
KonzoData.O.tr.status.f <- prune_taxa(f_0.0001, KonzoData.O.tr.status)

#write.csv((KonzoData.O@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Order_ReadCounts.csv")
#write.csv((KonzoData.O.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Order_RelAbund.csv")
#write.csv(t(KonzoData.O.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Order_Avg_RelAbund.csv")

#Mean and Standard Deviation
KonzoData.O.tr.df <- as.data.frame(t(KonzoData.O.tr@otu_table))
KonzoData.O.tr.df <- cbind(KonzoData.O.tr.df, KonzoData.O.tr@sam_data$Status)

colnames(KonzoData.O.tr.df)[colnames(KonzoData.O.tr.df)=="KonzoData.O.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.O.tr.df))
  {KonzoData.O.tr.df[i,]$Status <- KonzoData.O.tr@sam_data[rownames(KonzoData.O.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.O.tr.avg <- KonzoData.O.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.O.tr.avg.x <- t(KonzoData.O.tr.avg)
colnames(KonzoData.O.tr.avg.x) <- KonzoData.O.tr.avg.x[1,]   
KonzoData.O.tr.avg.x <- KonzoData.O.tr.avg.x[-1,]
colnames(KonzoData.O.tr.avg.x) <- paste("Avg", colnames(KonzoData.O.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.O.tr.sd <- KonzoData.O.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.O.tr.sd.x <- t(KonzoData.O.tr.sd)
colnames(KonzoData.O.tr.sd.x) <- KonzoData.O.tr.sd.x[1,]   
KonzoData.O.tr.sd.x <- KonzoData.O.tr.sd.x[-1,]
colnames(KonzoData.O.tr.sd.x) <- paste("SD", colnames(KonzoData.O.tr.sd.x), sep = "_")                                                 

KonzoData.O.tr.avg.sd <-merge(KonzoData.O.tr.avg.x,KonzoData.O.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.O.tr.avg.sd) <- KonzoData.O.tr.avg.sd[,1]   
KonzoData.O.tr.avg.sd <- KonzoData.O.tr.avg.sd[,-1]  
                           
KonzoData.O.tr.avg.sd <- KonzoData.O.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.O.tr.avg.sd, file = "./KonzoDataOrder_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.O.tr.avg.sd.f <- subset(KonzoData.O.tr.avg.sd, rownames(KonzoData.O.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.O.tr.avg.sd.f, file = "./KonzoDataOrder_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
Order.tr <- merge(KonzoData.O.tr.avg.sd,as.data.frame(KonzoData.O.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Order.tr, file = "./KonzoDataOrder_RelAbund_Supp.csv")                           

Order.tr.f <- merge(KonzoData.O.tr.avg.sd.f,as.data.frame(KonzoData.O.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Order.tr.f, file = "./KonzoDataOrder_RelAbund_filtered_Supp.csv")                           
                           

                           
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
KonzoData.F@sam_data$Status <- factor(KonzoData.F@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
#Read Counts to Relative Abundance
KonzoData.F.tr <- transform_sample_counts(KonzoData.F, function(x) x / sum(x))

KonzoData.F.tr.status <- merge_samples(KonzoData.F.tr, KonzoData.F.tr@sam_data$Status)
KonzoData.F.tr.status <- transform_sample_counts(KonzoData.F.tr.status, function(x) x / 30)                                          
#write.csv(t(KonzoData.F.tr.status@otu_table), file = "./KonzoDataFamily_AvgRelAbund_ByStatus.csv")
                                                
                                                 
#Filter
                                                                                                                        
Kinshasa.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa", KonzoData.F)
Kinshasa.F.tr <- transform_sample_counts(Kinshasa.F, function(x) x / sum(x))
Masimanimba.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba", KonzoData.F)
Masimanimba.F.tr <- transform_sample_counts(Masimanimba.F, function(x) x / sum(x))                                        
ULPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.F)
ULPZ.F.tr <- transform_sample_counts(ULPZ.F, function(x) x / sum(x))
KLPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.F)
KLPZ.F.tr <- transform_sample_counts(KLPZ.F, function(x) x / sum(x))
UHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.F)
UHPZ.F.tr <- transform_sample_counts(UHPZ.F, function(x) x / sum(x))
KHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)
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

#write.csv(filterList, file = "Kinshasa_Konzo3_Family_f_0.0001.csv")
                                 
x <- read.csv("Kinshasa_Konzo3_Family_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.F.f <- prune_taxa(f_0.0001, KonzoData.F)
KonzoData.F.tr.f <- prune_taxa(f_0.0001, KonzoData.F.tr)

KonzoData.F.tr.status.f <- prune_taxa(f_0.0001, KonzoData.F.tr.status)

#write.csv((KonzoData.F@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Family_ReadCounts.csv")
#write.csv((KonzoData.F.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Family_RelAbund.csv")
#write.csv(t(KonzoData.F.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Family_Avg_RelAbund.csv")

#Mean and Standard Deviation
KonzoData.F.tr.df <- as.data.frame(t(KonzoData.F.tr@otu_table))
KonzoData.F.tr.df <- cbind(KonzoData.F.tr.df, KonzoData.F.tr@sam_data$Status)

colnames(KonzoData.F.tr.df)[colnames(KonzoData.F.tr.df)=="KonzoData.F.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.F.tr.df))
  {KonzoData.F.tr.df[i,]$Status <- KonzoData.F.tr@sam_data[rownames(KonzoData.F.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.F.tr.avg <- KonzoData.F.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.F.tr.avg.x <- t(KonzoData.F.tr.avg)
colnames(KonzoData.F.tr.avg.x) <- KonzoData.F.tr.avg.x[1,]   
KonzoData.F.tr.avg.x <- KonzoData.F.tr.avg.x[-1,]
colnames(KonzoData.F.tr.avg.x) <- paste("Avg", colnames(KonzoData.F.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.F.tr.sd <- KonzoData.F.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.F.tr.sd.x <- t(KonzoData.F.tr.sd)
colnames(KonzoData.F.tr.sd.x) <- KonzoData.F.tr.sd.x[1,]   
KonzoData.F.tr.sd.x <- KonzoData.F.tr.sd.x[-1,]
colnames(KonzoData.F.tr.sd.x) <- paste("SD", colnames(KonzoData.F.tr.sd.x), sep = "_")                                                 

KonzoData.F.tr.avg.sd <-merge(KonzoData.F.tr.avg.x,KonzoData.F.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.F.tr.avg.sd) <- KonzoData.F.tr.avg.sd[,1]   
KonzoData.F.tr.avg.sd <- KonzoData.F.tr.avg.sd[,-1]  
                           
KonzoData.F.tr.avg.sd <- KonzoData.F.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.F.tr.avg.sd, file = "./KonzoDataFamily_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.F.tr.avg.sd.f <- subset(KonzoData.F.tr.avg.sd, rownames(KonzoData.F.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.F.tr.avg.sd.f, file = "./KonzoDataFamily_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
Family.tr <- merge(KonzoData.F.tr.avg.sd,as.data.frame(KonzoData.F.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Family.tr, file = "./KonzoDataFamily_RelAbund_Supp.csv")                           

Family.tr.f <- merge(KonzoData.F.tr.avg.sd.f,as.data.frame(KonzoData.F.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Family.tr.f, file = "./KonzoDataFamily_RelAbund_filtered_Supp.csv")                                
                           
                           
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
KonzoData.G@sam_data$Status <- factor(KonzoData.G@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
KonzoData.G@sam_data$Region <- factor(KonzoData.G@sam_data$Region, levels = c("Kinshasa", "Masimanimba", "Kahemba"))

#Read Count to Relative Abundance
KonzoData.G.tr <- transform_sample_counts(KonzoData.G, function(x) x / sum(x))
#Log Base 10 transform rel abund                                          
KonzoData.G.tr.log10 <- transform_sample_counts(KonzoData.G.tr, function(x) log10(x)) 

KonzoData.G.tr.status <- merge_samples(KonzoData.G.tr, KonzoData.G.tr@sam_data$Status, fun = mean)
KonzoData.G.tr.status <- transform_sample_counts(KonzoData.G.tr.status, function(x) x / 30)
 
                                               
#Filter                           
Kinshasa.G <- prune_samples(KonzoData.G@sam_data$Status == "Kinshasa", KonzoData.G)
Kinshasa.G.tr <- transform_sample_counts(Kinshasa.G, function(x) x / sum(x))
Masimanimba.G <- prune_samples(KonzoData.G@sam_data$Status == "Masimanimba", KonzoData.G)
Masimanimba.G.tr <- transform_sample_counts(Masimanimba.G, function(x) x / sum(x))                                        
ULPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.G)
ULPZ.G.tr <- transform_sample_counts(ULPZ.G, function(x) x / sum(x))
KLPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.G)
KLPZ.G.tr <- transform_sample_counts(KLPZ.G, function(x) x / sum(x))
UHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.G)
UHPZ.G.tr <- transform_sample_counts(UHPZ.G, function(x) x / sum(x))
KHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.G)
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


#write.csv(filterList, file = "Kinshasa_Konzo3_Genus_f_0.0001.csv")

x <- read.csv("Kinshasa_Konzo3_Genus_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.G.f <- prune_taxa(f_0.0001, KonzoData.G)
KonzoData.G.tr.f <- prune_taxa(f_0.0001, KonzoData.G.tr)

KonzoData.G.tr.status.f <- prune_taxa(f_0.0001, KonzoData.G.tr.status)

#write.csv((KonzoData.G@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Genus_ReadCounts.csv")
#write.csv((KonzoData.G.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Genus_RelAbund.csv")
#write.csv(t(KonzoData.G.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Genus_Avg_RelAbund.csv")

#Mean and Standard Deviation
KonzoData.G.tr.df <- as.data.frame(t(KonzoData.G.tr@otu_table))
KonzoData.G.tr.df <- cbind(KonzoData.G.tr.df, KonzoData.G.tr@sam_data$Status)

colnames(KonzoData.G.tr.df)[colnames(KonzoData.G.tr.df)=="KonzoData.G.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.G.tr.df))
  {KonzoData.G.tr.df[i,]$Status <- KonzoData.G.tr@sam_data[rownames(KonzoData.G.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.G.tr.avg <- KonzoData.G.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.G.tr.avg.x <- t(KonzoData.G.tr.avg)
colnames(KonzoData.G.tr.avg.x) <- KonzoData.G.tr.avg.x[1,]   
KonzoData.G.tr.avg.x <- KonzoData.G.tr.avg.x[-1,]
colnames(KonzoData.G.tr.avg.x) <- paste("Avg", colnames(KonzoData.G.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.G.tr.sd <- KonzoData.G.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.G.tr.sd.x <- t(KonzoData.G.tr.sd)
colnames(KonzoData.G.tr.sd.x) <- KonzoData.G.tr.sd.x[1,]   
KonzoData.G.tr.sd.x <- KonzoData.G.tr.sd.x[-1,]
colnames(KonzoData.G.tr.sd.x) <- paste("SD", colnames(KonzoData.G.tr.sd.x), sep = "_")                                                 

KonzoData.G.tr.avg.sd <-merge(KonzoData.G.tr.avg.x,KonzoData.G.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.G.tr.avg.sd) <- KonzoData.G.tr.avg.sd[,1]   
KonzoData.G.tr.avg.sd <- KonzoData.G.tr.avg.sd[,-1]  
                           
KonzoData.G.tr.avg.sd <- KonzoData.G.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.G.tr.avg.sd, file = "./KonzoDataGenus_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.G.tr.avg.sd.f <- subset(KonzoData.G.tr.avg.sd, rownames(KonzoData.G.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.G.tr.avg.sd.f, file = "./KonzoDataGenus_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
Genus.tr <- merge(KonzoData.G.tr.avg.sd,as.data.frame(KonzoData.G.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Genus.tr, file = "./KonzoDataGenus_RelAbund_Supp.csv")                           

Genus.tr.f <- merge(KonzoData.G.tr.avg.sd.f,as.data.frame(KonzoData.G.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Genus.tr.f, file = "./KonzoDataGenus_RelAbund_filtered_Supp.csv")                                
                           
                                                   
#Needed later for Figure 3 (Geography excluding all konzo individuals)                           
Geography.G <- prune_samples((KonzoData.G@sam_data$Status != "Konzo_Low_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_High_Prevalence_Zone"), KonzoData.G)                                              
Geography.G.tr <-  transform_sample_counts(Geography.G, function(x) x / sum(x))
Geography.G.tr.log10 <-  transform_sample_counts(Geography.G.tr, function(x) log10(x))    
Geography.G.f <- prune_taxa(f_0.0001, Geography.G)                                                 
Geography.G.tr.f <- prune_taxa(f_0.0001, Geography.G.tr)
                           
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

KonzoData.S@sam_data$Status <- factor(KonzoData.S@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

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
ULPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.S)
ULPZ.S.tr <- transform_sample_counts(ULPZ.S, function(x) x / sum(x))
KLPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.S)
KLPZ.S.tr <- transform_sample_counts(KLPZ.S, function(x) x / sum(x))
UHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.S)
UHPZ.S.tr <- transform_sample_counts(UHPZ.S, function(x) x / sum(x))
KHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)
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

#write.csv(filterList, file = "Kinshasa_Konzo3_Species_f_0.0001.csv")
#Manually Change Streptococcus sp. 'group B' to Streptococcus sp. group B 
#Otherwise the filter function thinks its a different species and will filter it out                           

x <- read.csv("Kinshasa_Konzo3_Species_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.S.f <- prune_taxa(f_0.0001, KonzoData.S)
KonzoData.S.tr.f <- prune_taxa(f_0.0001, KonzoData.S.tr)

KonzoData.S.tr.status.f <- prune_taxa(f_0.0001, KonzoData.S.tr.status)
                                                                                                       
#write.csv((KonzoData.S.tr@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Species_RelAbund.csv")
#write.csv((KonzoData.S@otu_table), file = "./KonzoMicrobiome_Samples_Bacteria_Species_ReadCounts.csv")
#write.csv(t(KonzoData.S.tr.status@otu_table), file = "./KonzoMicrobiome_Groups_Bacteria_Species_Avg_RelAbund.csv")
                                                 
#Mean and Standard Deviation
KonzoData.S.tr.df <- as.data.frame(t(KonzoData.S.tr@otu_table))
KonzoData.S.tr.df <- cbind(KonzoData.S.tr.df, KonzoData.S.tr@sam_data$Status)

colnames(KonzoData.S.tr.df)[colnames(KonzoData.S.tr.df)=="KonzoData.S.tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData.S.tr.df))
  {KonzoData.S.tr.df[i,]$Status <- KonzoData.S.tr@sam_data[rownames(KonzoData.S.tr.df[i,]),]$Status
  } 
                                                 
KonzoData.S.tr.avg <- KonzoData.S.tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData.S.tr.avg.x <- t(KonzoData.S.tr.avg)
colnames(KonzoData.S.tr.avg.x) <- KonzoData.S.tr.avg.x[1,]   
KonzoData.S.tr.avg.x <- KonzoData.S.tr.avg.x[-1,]
colnames(KonzoData.S.tr.avg.x) <- paste("Avg", colnames(KonzoData.S.tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData.S.tr.sd <- KonzoData.S.tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData.S.tr.sd.x <- t(KonzoData.S.tr.sd)
colnames(KonzoData.S.tr.sd.x) <- KonzoData.S.tr.sd.x[1,]   
KonzoData.S.tr.sd.x <- KonzoData.S.tr.sd.x[-1,]
colnames(KonzoData.S.tr.sd.x) <- paste("SD", colnames(KonzoData.S.tr.sd.x), sep = "_")                                                 

KonzoData.S.tr.avg.sd <-merge(KonzoData.S.tr.avg.x,KonzoData.S.tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData.S.tr.avg.sd) <- KonzoData.S.tr.avg.sd[,1]   
KonzoData.S.tr.avg.sd <- KonzoData.S.tr.avg.sd[,-1]  
                           
KonzoData.S.tr.avg.sd <- KonzoData.S.tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData.S.tr.avg.sd, file = "./KonzoDataSpecies_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData.S.tr.avg.sd.f <- subset(KonzoData.S.tr.avg.sd, rownames(KonzoData.S.tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData.S.tr.avg.sd.f, file = "./KonzoDataSpecies_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
Species.tr <- merge(KonzoData.S.tr.avg.sd,as.data.frame(KonzoData.S.tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(Species.tr, file = "./KonzoDataSpecies_RelAbund_Supp.csv")                           

Species.tr.f <- merge(KonzoData.S.tr.avg.sd.f,as.data.frame(KonzoData.S.tr.f@otu_table),by='row.names', sort = FALSE)                           
write.csv(Species.tr.f, file = "./KonzoDataSpecies_RelAbund_filtered_Supp.csv")
                           
### Estimate Richness
                           
#for any sample where the relabund is <0.0001, set the read count to 0 for alpha diversity calcs                           
setzero = function(x){
  x[(x / sum(x)) < (1e-4)] <- 0
  return(x)
}                           
#Transform Read Count from KonzoData.S (Bacteria Species data)                           
KonzoData.S.0 <- transform_sample_counts(KonzoData.S, setzero)
                           
otuD.S.0 <- as.data.frame(t(otu_table(KonzoData.S.0)))
diversity.S.0 <- estimate_richness(KonzoData.S.0)
diversity.S.0 <- cbind(sample_data(KonzoData.S.0),diversity.S.0) #Check if correct sample data was cbind. Can be tricky so always confirm
diversity.S.0$Status <- as.factor(diversity.S.0$Status)
diversity.S.0$Status <- factor(diversity.S.0$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

#Shapiro-Wilk Normality Test
shapiro.test(diversity.S.0$Observed) #p-value = 0.09588
#shapiro.test(diversity.S.0$Chao1)                           
shapiro.test(diversity.S.0$Shannon) # p-value = 1.345e-09                           
#shapiro.test(diversity.S.0$ACE)
shapiro.test(diversity.S.0$Simpson) #p-value < 2.2e-16
shapiro.test(diversity.S.0$Fisher) #p-value = 0.3341
                           
#methods that are reliant on singletons such as ACE cannot be accurately calculated. Here only methods not reliant on singletons are assessed, so warning can be ignored                           
#STATISTICS for Estimate Richness
#One-way ANOVA (or Kruskal) to see if there is a statitically significant difference in the measure of alpha diversity and output saved in txt file
observed.aov <- aov(Observed ~ Status, data = diversity.S.0) #p-value = 0.00717
fisher.aov <- aov(Fisher ~ Status, data = diversity.S.0) #p-value = 0.00762
                           
shannon.kru <- kruskal.test(Shannon ~ Status, data = diversity.S.0) #p-value = 0.02398515
simpson.kru <- kruskal.test(Simpson ~ Status, data = diversity.S.0)#p-value = 0.01411952

observed.tukey <- TukeyHSD(observed.aov, data = diversity.S.0)
fisher.tukey <- TukeyHSD(fisher.aov, data = diversity.S.0)
 
shannon.dunn <- dunnTest(Shannon ~Status, data = diversity.S.0, method = "bh")  
simpson.dunn <- dunnTest(Simpson ~Status, data = diversity.S.0, method = "bh")                             
                           
shannon.wilcox <- pairwise.wilcox.test(diversity.S.0$Shannon, diversity.S.0$Status, data = diversity.S.0, p.adjust.method = "BH")
simpson.wilcox <- pairwise.wilcox.test(diversity.S.0$Simpson, diversity.S.0$Status, data = diversity.S.0, p.adjust.method = "BH")                          

                                                     
write("Observed ~ Status ANOVA test", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(summary(observed.aov), append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("TukeyHSA observed.aov", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(observed.tukey, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 
                           
write("Shannon ~ Status Kruskal test", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(shannon.kru, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("Post-hoc shannon.dunn", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(shannon.dunn, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("Pairewise MWW shannon.kru", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(shannon.wilcox, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("Simpson ~ Status Kruskal test", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(simpson.kru, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 
                           
write("Post-hoc simpson.dunn", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(simpson.dunn, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("Pairewise MWW simpson.kru", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(simpson.wilcox, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt")                            

write("Fisher ~ Status ANOVA test", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(summary(fisher.aov), append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt") 

write("TukeyHSA fisher.aov", file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt" ,append=TRUE)
capture.output(fisher.tukey, append = TRUE, file="KinshasaControl_Konzo3_Bacteria_Species_SetZeroData_EstimateRichness_Statistics.txt")
                           
                           
###Figure 2 ------------------------------ 
                           
observed <- ggplot(diversity.S.0, aes(factor(Status), Observed)) + geom_boxplot(aes(fill = factor(Status)),fatten = 1, outlier.shape = NA) + labs(x = element_blank(), y = "Species") + theme(axis.text.x = element_blank()) + theme_classic()
observed <- observed + geom_jitter(position=position_jitter(0.2), size = 0.3)
observed2 <- observed + stat_summary(fun=mean, geom="point", shape=23, size=1.5, color = "black", fill="white")
observed3 <- observed2 + theme(legend.position="bottom", legend.margin=margin(-10,0,0,0)) + theme(legend.direction = "horizontal") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(axis.ticks.x = element_blank(), axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_blank())

observed4 <- observed3 + guides(fill=guide_legend(ncol=6)) 
observed4 <- observed4 + scale_fill_manual(values = konzo_color, labels = SSSL)

shan <- ggplot(diversity.S.0, aes(factor(Status), Shannon))+ geom_boxplot(aes(fill = factor(Status)),fatten = 1, outlier.shape = NA) + labs(x = element_blank(), y = "Shannon Diversity Index") + theme(axis.text.x = element_blank()) + theme_classic()
shan <- shan + geom_jitter(position=position_jitter(0.2), size = 0.3)
shan2 <- shan + stat_summary(fun=mean, geom="point", shape=23, size=1.5, color = "black", fill="white")
shan3 <- shan2 + theme(legend.position="bottom", legend.margin=margin(-10,0,0,0)) + theme(legend.direction = "horizontal") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 7), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(axis.ticks.x = element_blank(), axis.title.y = element_text(size = 5), axis.text.y = element_text(size = 7), axis.text.x = element_blank())

shan4 <- shan3 + guides(fill=guide_legend(ncol=6)) 
shan4 <- shan4 + scale_fill_manual(values = konzo_color, labels = SSSL)

#To get Top taxa
                                                    
                           
#Phylum
                           
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Phylum")
                           
Kinshasa.P <- prune_samples(KonzoData.P@sam_data$Status == "Kinshasa", KonzoData.P)
Kinshasa.P.tr <- transform_sample_counts(Kinshasa.P, function(x) x / sum(x))
Masimanimba.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba", KonzoData.P)
Masimanimba.P.tr <- transform_sample_counts(Masimanimba.P, function(x) x / sum(x))                                        
ULPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.P)
ULPZ.P.tr <- transform_sample_counts(ULPZ.P, function(x) x / sum(x))
KLPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.P)
KLPZ.P.tr <- transform_sample_counts(KLPZ.P, function(x) x / sum(x))
UHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.P)
UHPZ.P.tr <- transform_sample_counts(UHPZ.P, function(x) x / sum(x))
KHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
KHPZ.P.tr <- transform_sample_counts(KHPZ.P, function(x) x / sum(x))
   

top1 = sort(taxa_sums(Kinshasa.P.tr), TRUE)[1:4]
top2 = sort(taxa_sums(Masimanimba.P.tr), TRUE)[1:4]
top3 = sort(taxa_sums(ULPZ.P.tr), TRUE)[1:4]
top4 = sort(taxa_sums(KLPZ.P.tr), TRUE)[1:4]
top5 = sort(taxa_sums(UHPZ.P.tr), TRUE)[1:4]
top6 = sort(taxa_sums(KHPZ.P.tr), TRUE)[1:4]
                                                                          
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_P <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ
                                     
write.csv(top_P, file = "Kinshasa_Konzo3_Phylum_Top4.csv")

#Class 
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Class")
                                     
Kinshasa.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa", KonzoData.C)
Kinshasa.C.tr <- transform_sample_counts(Kinshasa.C, function(x) x / sum(x))
Masimanimba.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba", KonzoData.C)
Masimanimba.C.tr <- transform_sample_counts(Masimanimba.C, function(x) x / sum(x))                                        
ULPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.C)
ULPZ.C.tr <- transform_sample_counts(ULPZ.C, function(x) x / sum(x))
KLPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.C)
KLPZ.C.tr <- transform_sample_counts(KLPZ.C, function(x) x / sum(x))
UHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.C)
UHPZ.C.tr <- transform_sample_counts(UHPZ.C, function(x) x / sum(x))
KHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)
KHPZ.C.tr <- transform_sample_counts(KHPZ.C, function(x) x / sum(x))

top1 = sort(taxa_sums(Kinshasa.C.tr), TRUE)[1:5]
top2 = sort(taxa_sums(Masimanimba.C.tr), TRUE)[1:5]
top3 = sort(taxa_sums(ULPZ.C.tr), TRUE)[1:5]
top4 = sort(taxa_sums(KLPZ.C.tr), TRUE)[1:5]
top5 = sort(taxa_sums(UHPZ.C.tr), TRUE)[1:5]
top6 = sort(taxa_sums(KHPZ.C.tr), TRUE)[1:5]

                                                                          
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_C <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(top_C, file = "Kinshasa_Konzo3_Class_Top5.csv")                                     
                                     
#ORDER
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Order")
                                     
Kinshasa.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa", KonzoData.O)
Kinshasa.O.tr <- transform_sample_counts(Kinshasa.O, function(x) x / sum(x))
Masimanimba.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba", KonzoData.O)
Masimanimba.O.tr <- transform_sample_counts(Masimanimba.O, function(x) x / sum(x))                                        
ULPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.O)
ULPZ.O.tr <- transform_sample_counts(ULPZ.O, function(x) x / sum(x))
KLPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_Konzo_Prevalence_Zone", KonzoData.O)
KLPZ.O.tr <- transform_sample_counts(KLPZ.O, function(x) x / sum(x))
UHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.O)
UHPZ.O.tr <- transform_sample_counts(UHPZ.O, function(x) x / sum(x))
KHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)
KHPZ.O.tr <- transform_sample_counts(KHPZ.O, function(x) x / sum(x))

top1 = sort(taxa_sums(Kinshasa.O.tr), TRUE)[1:5]
top2 = sort(taxa_sums(Masimanimba.O.tr), TRUE)[1:5]
top3 = sort(taxa_sums(ULPZ.O.tr), TRUE)[1:5]
top4 = sort(taxa_sums(KLPZ.O.tr), TRUE)[1:5]
top5 = sort(taxa_sums(UHPZ.O.tr), TRUE)[1:5]
top6 = sort(taxa_sums(KHPZ.O.tr), TRUE)[1:5]
                                     
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_O <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ    
                                     
write.csv(top_O, file = "Kinshasa_Konzo3_Order_Top5.csv")
                                     
#FAMILY
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Family")
                           
Kinshasa.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa", KonzoData.F)
Kinshasa.F.tr <- transform_sample_counts(Kinshasa.F, function(x) x / sum(x))
Masimanimba.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba", KonzoData.F)
Masimanimba.F.tr <- transform_sample_counts(Masimanimba.F, function(x) x / sum(x))                                        
ULPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.F)
ULPZ.F.tr <- transform_sample_counts(ULPZ.F, function(x) x / sum(x))
KLPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.F)
KLPZ.F.tr <- transform_sample_counts(KLPZ.F, function(x) x / sum(x))
UHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.F)
UHPZ.F.tr <- transform_sample_counts(UHPZ.F, function(x) x / sum(x))
KHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)
KHPZ.F.tr <- transform_sample_counts(KHPZ.F, function(x) x / sum(x))
                          
top1 = sort(taxa_sums(Kinshasa.F.tr), TRUE)[1:5]
top2 = sort(taxa_sums(Masimanimba.F.tr), TRUE)[1:5]
top3 = sort(taxa_sums(ULPZ.F.tr), TRUE)[1:5]
top4 = sort(taxa_sums(KLPZ.F.tr), TRUE)[1:5]
top5 = sort(taxa_sums(UHPZ.F.tr), TRUE)[1:5]
top6 = sort(taxa_sums(KHPZ.F.tr), TRUE)[1:5]
                                     
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_F <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(top_F, file = "Kinshasa_Konzo3_Family_Top5.csv")
                                     
#GENUS
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")
                                     
Kinshasa.G <- prune_samples(KonzoData.G@sam_data$Status == "Kinshasa", KonzoData.G)
Kinshasa.G.tr <- transform_sample_counts(Kinshasa.G, function(x) x / sum(x))
Masimanimba.G <- prune_samples(KonzoData.G@sam_data$Status == "Masimanimba", KonzoData.G)
Masimanimba.G.tr <- transform_sample_counts(Masimanimba.G, function(x) x / sum(x))                                        
ULPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.G)
ULPZ.G.tr <- transform_sample_counts(ULPZ.G, function(x) x / sum(x))
KLPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.G)
KLPZ.G.tr <- transform_sample_counts(KLPZ.G, function(x) x / sum(x))
UHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.G)
UHPZ.G.tr <- transform_sample_counts(UHPZ.G, function(x) x / sum(x))
KHPZ.G <- prune_samples(KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.G)
KHPZ.G.tr <- transform_sample_counts(KHPZ.G, function(x) x / sum(x))                           

top1 = sort(taxa_sums(Kinshasa.G.tr), TRUE)[1:7]
top2 = sort(taxa_sums(Masimanimba.G.tr), TRUE)[1:7]
top3 = sort(taxa_sums(ULPZ.G.tr), TRUE)[1:7]
top4 = sort(taxa_sums(KLPZ.G.tr), TRUE)[1:7]
top5 = sort(taxa_sums(UHPZ.G.tr), TRUE)[1:7]
top6 = sort(taxa_sums(KHPZ.G.tr), TRUE)[1:7]
                                     
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_G <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(top_G, file = "Kinshasa_Konzo3_Genus_Top7.csv")
                                     
#SPECIES
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Species")
                                     
Kinshasa.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa", KonzoData.S)
Kinshasa.S.tr <- transform_sample_counts(Kinshasa.S, function(x) x / sum(x))
Masimanimba.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba", KonzoData.S)
Masimanimba.S.tr <- transform_sample_counts(Masimanimba.S, function(x) x / sum(x))                                        
ULPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.S)
ULPZ.S.tr <- transform_sample_counts(ULPZ.S, function(x) x / sum(x))
KLPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.S)
KLPZ.S.tr <- transform_sample_counts(KLPZ.S, function(x) x / sum(x))
UHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.S)
UHPZ.S.tr <- transform_sample_counts(UHPZ.S, function(x) x / sum(x))
KHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)
KHPZ.S.tr <- transform_sample_counts(KHPZ.S, function(x) x / sum(x))

top1 = sort(taxa_sums(Kinshasa.S.tr), TRUE)[1:20]
top2 = sort(taxa_sums(Masimanimba.S.tr), TRUE)[1:20]
top3 = sort(taxa_sums(ULPZ.S.tr), TRUE)[1:20]
top4 = sort(taxa_sums(KLPZ.S.tr), TRUE)[1:20]
top5 = sort(taxa_sums(UHPZ.S.tr), TRUE)[1:20]
top6 = sort(taxa_sums(KHPZ.S.tr), TRUE)[1:20]
                                     
top12 <- union(names(top1),names(top2)) #Kin, Mas
top34 <- union(names(top3), names(top4)) #ULPZ, KLPZ
top56 <- union(names(top5), names(top6))
top1234 <- union(top12, top34) #Kin, Mas, ULPZ, KLPZ
top_S <- union(top1234, top56) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

write.csv(top_S, file = "Kinshasa_Konzo3_Species_Top20.csv")                           
                           
#Stacked Bar plots
#Genus
#Perhaps a rearrangement of the colorscheme for this                                     
#t3_cols = c("#00B9E3","#D39200","#00C19F","#F8766D","#619CFF","#93AA00", "#DB72FB", "#00BA38", "#FF61C3")
                                     
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")
                                     
top_G <- read.csv("Kinshasa_Konzo3_Genus_Top7.csv", row.names = 1, colClasses = "character")
top_G <- unlist(top_G)

KonzoData.G.tr.top = prune_taxa(top_G, KonzoData.G.tr)
physeqdf <- psmelt(KonzoData.G.tr.top)
p <- ggplot(physeqdf, aes(x=Abundance, y=reorder(Sample, -Abundance), fill = reorder(genus, Abundance)))
p <- p + geom_bar(stat="identity", width = 1)  
p$data$Status <- factor(p$data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))                               
#p$data$genus <- factor(p$data$genus, levels = c("Prevotella","Bacteroides","Faecalibacterium","Escherichia","Alistipes", "Eubacterium", "Bifidobacterium","Butyricimonas","Anaerostipes") )
p <- p + labs(y = element_blank(), x = "Relative Abundance") +  scale_fill_discrete(name = "Genus")
top_genus_plot <- p + theme(legend.position="rigth") + theme(legend.key.size = unit(.2, "cm"), legend.text = element_text(face="italic"))

top_genus_plot <- top_genus_plot + facet_grid(rows = vars(Status), scales = "free_y", space = "free", switch = "y", labeller = as_labeller(SSSL)) +
  theme(strip.background = element_blank(), strip.placement = "outside", strip.text.y.left = element_text(angle = 0, vjust=0.5, hjust=0, size = 7))

top_genus_plot<- top_genus_plot + theme(panel.spacing=unit(0, "lines"), panel.border = element_rect(color = "black", fill = NA, size = 0.5))

top_genus_plot <- top_genus_plot + 
  theme(legend.position="bottom", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,0,0)) + 
  scale_x_continuous(expand = c(0,0), limits = c(0,0.9)) +
  theme(plot.title = element_blank(), legend.key.size = unit(.2, "cm"), legend.text = element_text(size = 5), legend.title = element_text(size = 5)) + 
  guides(fill=guide_legend(ncol=3,byrow=TRUE)) + # scale_fill_manual(values = t3_cols) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x = element_text(size = 7), axis.title.x = element_text(size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
top_genus_plot

ad <- ggarrange(observed4, shan4, labels = c("A", "B"), font.label = list(size = 7), ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", align = "hv") 

Gen <-  ggarrange(top_genus_plot, labels = c("C"), font.label = list(size = 7), ncol = 1, nrow = 1) 

s <- plot_spacer() + theme_minimal()

placeholder <-  ggarrange(s, labels = c("D"), font.label = list(size = 7), ncol = 1, nrow = 1)    

Gen_ph <- ggarrange(Gen,placeholder, widths = c(1, 1), ncol = 2, nrow = 1)

F1 <- arrangeGrob(ad, Gen_ph, ncol = 1, nrow = 2,
             layout_matrix = rbind(c(1), c(2), c(2), c(2)))
#Figure 2
tiff(filename = "KinshasaKonzo3_TaxaFigure_WithoutHeatMap.tiff", width = 5, height = 5, units = "in", res = 600)
as_ggplot(F1)
dev.off()
                                     
                                     
#When adding heat map in gimp to full figure, make sure Image > Print Size has correct inches and ppi (set to >=300))
o <- as.data.frame(otu_table(KonzoData.S.tr.status.f))                                                 
tiff(filename = "KinshasaKonzo3_Bacteria_Species_Heatmap_V1.tiff", width = 3.5, height = 2.8, units = "in", res = 600)
heatmap.2(as.matrix(t(o)), scale = "row", trace = "none", keysize = 0.25, labRow = "Species", labCol = SSSL, margins = c(1, 1), Rowv = FALSE, dendrogram = "column", key.title = NA, srtCol = 0, srtRow = 90 , cexCol = 0.75, cexRow = 0.75, offsetRow = 0, offsetCol = 0, lhei = c(0.5,2,2,2), lwid = c(0.1,1,1,1), key.par = list(cex=0.5), lmat = rbind(c(0,3,3,0),c(2,1,1,0),c(2,1,1,0),c(2,1,1,4)), adjCol = c(0.5,0.5), adjRow = c(4.5,0.25))
dev.off() 
                                     
#When adding heat map in gimp to full figure, make sure Image > Print Size has correct inches and ppi (set to >=300))
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Species")

o <- as.data.frame(otu_table(KonzoData.S.tr.status.f))                                                 
tiff(filename = "KinshasaKonzo3_Bacteria_Species_Heatmap_V2.tiff", width = 2.25, height = 3.75, units = "in", res = 600)
heatmap.2(as.matrix(t(o)), scale = "row", trace = "none", keysize = 0.25, labRow = "Species", labCol = SSSL, margins = c(1, 1), Rowv = FALSE, dendrogram = "column", key.title = NA, srtCol = 0, srtRow = 90 , cexCol = 0.75, cexRow = 0.75, offsetRow = 0, offsetCol = 0, lhei = c(0.5,2,2,1.25), lwid = c(0.1,1,1), key.par = list(cex=0.5), lmat = rbind(c(0,3,3),c(2,1,1),c(2,1,1),c(0,0,4)), adjCol = c(0.5,0.5), adjRow = c(4.5,0.25))
dev.off()                                        
                                                      
#------------------------------------------------
     
### ALDEx2 : aldex is used to convert raw read counts (aftering filtering out low abundance taxa) to CLR to infer abuandance and perform statistical testing. 
#Aldex does both he Welch's t-test and Wilcoxon test on the pairewise comparisions and doing BH correction on the expected values, but only Wilcoxon test results are reported in figures etc.                                   

#Supplemental File 3 where wi.ep and wi.eBH are reported at each taxa (the final AD13 dataframe is in each sheet for the specific taxa rank)
                                     
#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

#PHYLUM
x <- read.csv("Kinshasa_Konzo3_Phylum_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

KonzoData.P.f <- prune_taxa(f_0.0001, KonzoData.P) #filtered readcount phyloseq object
                                               
#Urban vs Rural
P.aldex.DF <- as.data.frame(KonzoData.P.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.P.Aldex.DF <- x
KonzoData.P.CLR.DF <- t(x)
KonzoData.P.CLR.DF <- KonzoData.P.CLR.DF[4:183,]
rownames(KonzoData.P.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.P.CLR.DF))
KonzoData.P.CLR.DF <- cbind(KonzoData.P.CLR.DF, as.data.frame(KonzoData.P.f@sam_data$Status))
colnames(KonzoData.P.CLR.DF)[colnames(KonzoData.P.CLR.DF)=="KonzoData.P.f@sam_data$Status"] <- "Status"
write.csv(KonzoData.P.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Phylum_KonzoData_Aldex_MedianCLRValues.csv")  

conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  

#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.P <-  prune_samples((KonzoData.P@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.P@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.P)
KinMasLPZ.P.f <- prune_taxa(f_0.0001, KinMasLPZ.P)   

P.aldex.DF <- as.data.frame(KinMasLPZ.P.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.P <-  prune_samples((KonzoData.P@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.P@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.P)
KinMasHPZ.P.f <- prune_taxa(f_0.0001, KinMasHPZ.P)   

P.aldex.DF <- as.data.frame(KinMasHPZ.P.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_HPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.P <-  prune_samples((KonzoData.P@sam_data$Status == "Kinshasa") | (KonzoData.P@sam_data$Status == "Masimanimba"),  KonzoData.P)
KinMas.P.f <- prune_taxa(f_0.0001, KinMas.P)   

P.aldex.DF <- as.data.frame(KinMas.P.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Kin_vs_Mas_Aldex.csv")  



#Kin vs ULPZ
KinULPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Kinshasa") | (KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.P)
KinULPZ.P.f <- prune_taxa(f_0.0001, KinULPZ.P)   

P.aldex.DF <- as.data.frame(KinULPZ.P.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs ULPZ
MasULPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Masimanimba") | (KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.P)
MasULPZ.P.f <- prune_taxa(f_0.0001, MasULPZ.P)   

P.aldex.DF <- as.data.frame(MasULPZ.P.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs UHPZ
KinUHPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Kinshasa") | (KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.P)
KinUHPZ.P.f <- prune_taxa(f_0.0001, KinUHPZ.P)   

P.aldex.DF <- as.data.frame(KinUHPZ.P.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs UHPZ
MasUHPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Masimanimba") | (KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.P)
MasUHPZ.P.f <- prune_taxa(f_0.0001, MasUHPZ.P)   

P.aldex.DF <- as.data.frame(MasUHPZ.P.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Kinshasa") | (KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.P)
KinKLPZ.P.f <- prune_taxa(f_0.0001, KinKLPZ.P)   

P.aldex.DF <- as.data.frame(KinKLPZ.P.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Kin_vs_KLPZ_Aldex.csv")  

#Mas vs KLPZ
MasKLPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Masimanimba") | (KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.P)
MasKLPZ.P.f <- prune_taxa(f_0.0001, MasKLPZ.P)   

P.aldex.DF <- as.data.frame(MasKLPZ.P.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs KHPZ
KinKHPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Kinshasa") | (KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.P)
KinKHPZ.P.f <- prune_taxa(f_0.0001, KinKHPZ.P)   

P.aldex.DF <- as.data.frame(KinKHPZ.P.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs KHPZ
MasKHPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Masimanimba") | (KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.P)
MasKHPZ.P.f <- prune_taxa(f_0.0001, MasKHPZ.P)   

P.aldex.DF <- as.data.frame(MasKHPZ.P.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_Mas_vs_KHPZ_Aldex.csv")  

#ULPZ vs UHPZ

Control.P <-  prune_samples((KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.P)
Control.P.f <- prune_taxa(f_0.0001, Control.P)   

P.aldex.DF <- as.data.frame(Control.P.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_ULPZ_vs_UHPZ_Aldex.csv")  


#KLPZ vs. KHPZ
Disease.P <-  prune_samples((KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.P)
Disease.P.f <- prune_taxa(f_0.0001, Disease.P)   

P.aldex.DF <- as.data.frame(Disease.P.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30)) #Always Check this
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_KLPZ_vs_KHPZ_Aldex.csv")  

#ULPZ vs. KLPZ

LPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.P)
LPZ.P.f <- prune_taxa(f_0.0001, LPZ.P)   

P.aldex.DF <- as.data.frame(LPZ.P.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_ULPZ_vs_KLPZ_Aldex.csv")  

#UHPZ vs. KHPZ

HPZ.P <-  prune_samples((KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.P)
HPZ.P.f <- prune_taxa(f_0.0001, HPZ.P)   

P.aldex.DF <- as.data.frame(HPZ.P.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30)) #Always Check this
x <- aldex(P.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Phylum_UHPZ_vs_KHPZ_Aldex.csv")  
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Phylum_Filtered_Aldex.csv")  

                                     
##CLASS
x <- read.csv("Kinshasa_Konzo3_Class_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KonzoData
C.aldex.DF <- as.data.frame(KonzoData.C.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.C.Aldex.DF <- x
KonzoData.C.CLR.DF <- t(x)
KonzoData.C.CLR.DF <- KonzoData.C.CLR.DF[4:183,]
rownames(KonzoData.C.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.C.CLR.DF))
KonzoData.C.CLR.DF <- cbind(KonzoData.C.CLR.DF, as.data.frame(KonzoData.C.f@sam_data$Status))
colnames(KonzoData.C.CLR.DF)[colnames(KonzoData.C.CLR.DF)=="KonzoData.C.f@sam_data$Status"] <- "Status"
write.csv(KonzoData.C.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Class_KonzoData_Aldex_MedianCLRValues.csv")  

conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  


#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.C <-  prune_samples((KonzoData.C@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.C@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.C)
KinMasLPZ.C.f <- prune_taxa(f_0.0001, KinMasLPZ.C)   

C.aldex.DF <- as.data.frame(KinMasLPZ.C.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.C <-  prune_samples((KonzoData.C@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.C@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.C)
KinMasHPZ.C.f <- prune_taxa(f_0.0001, KinMasHPZ.C)   

C.aldex.DF <- as.data.frame(KinMasHPZ.C.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_HPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.C <-  prune_samples((KonzoData.C@sam_data$Status == "Kinshasa") | (KonzoData.C@sam_data$Status == "Masimanimba"),  KonzoData.C)
KinMas.C.f <- prune_taxa(f_0.0001, KinMas.C)   

C.aldex.DF <- as.data.frame(KinMas.C.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Kin_vs_Mas_Aldex.csv")  

#Kin vs ULPZ
KinULPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Kinshasa") | (KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.C)
KinULPZ.C.f <- prune_taxa(f_0.0001, KinULPZ.C)   

C.aldex.DF <- as.data.frame(KinULPZ.C.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs ULPZ
MasULPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Masimanimba") | (KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.C)
MasULPZ.C.f <- prune_taxa(f_0.0001, MasULPZ.C)   

C.aldex.DF <- as.data.frame(MasULPZ.C.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs UHPZ
KinUHPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Kinshasa") | (KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.C)
KinUHPZ.C.f <- prune_taxa(f_0.0001, KinUHPZ.C)   

C.aldex.DF <- as.data.frame(KinUHPZ.C.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs UHPZ
MasUHPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Masimanimba") | (KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.C)
MasUHPZ.C.f <- prune_taxa(f_0.0001, MasUHPZ.C)   

C.aldex.DF <- as.data.frame(MasUHPZ.C.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Kinshasa") | (KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.C)
KinKLPZ.C.f <- prune_taxa(f_0.0001, KinKLPZ.C)   

C.aldex.DF <- as.data.frame(KinKLPZ.C.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Kin_vs_KLPZ_Aldex.csv")  

#Mas vs KLPZ
MasKLPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Masimanimba") | (KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.C)
MasKLPZ.C.f <- prune_taxa(f_0.0001, MasKLPZ.C)   

C.aldex.DF <- as.data.frame(MasKLPZ.C.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs KHPZ
KinKHPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Kinshasa") | (KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.C)
KinKHPZ.C.f <- prune_taxa(f_0.0001, KinKHPZ.C)   

C.aldex.DF <- as.data.frame(KinKHPZ.C.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs KHPZ
MasKHPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Masimanimba") | (KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.C)
MasKHPZ.C.f <- prune_taxa(f_0.0001, MasKHPZ.C)   

C.aldex.DF <- as.data.frame(MasKHPZ.C.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_Mas_vs_KHPZ_Aldex.csv")  

#ULPZ vs UHPZ

Control.C <-  prune_samples((KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.C)
Control.C.f <- prune_taxa(f_0.0001, Control.C)   

C.aldex.DF <- as.data.frame(Control.C.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_ULPZ_vs_UHPZ_Aldex.csv")  

#KLPZ vs. KHPZ
Disease.C <-  prune_samples((KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.C)
Disease.C.f <- prune_taxa(f_0.0001, Disease.C)   

C.aldex.DF <- as.data.frame(Disease.C.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30)) #Always Check this
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_KLPZ_vs_KHPZ_Aldex.csv")  

#ULPZ vs. KLPZ

LPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.C)
LPZ.C.f <- prune_taxa(f_0.0001, LPZ.C)   

C.aldex.DF <- as.data.frame(LPZ.C.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_ULPZ_vs_KLPZ_Aldex.csv")  

#UHPZ vs. KHPZ

HPZ.C <-  prune_samples((KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.C)
HPZ.C.f <- prune_taxa(f_0.0001, HPZ.C)   

C.aldex.DF <- as.data.frame(HPZ.C.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30)) #Always Check this
x <- aldex(C.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Class_UHPZ_vs_KHPZ_Aldex.csv")  
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Class_Filtered_Aldex.csv")  

##ORDER
x <- read.csv("Kinshasa_Konzo3_Order_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KonzoData
O.aldex.DF <- as.data.frame(KonzoData.O.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.O.Aldex.DF <- x
KonzoData.O.CLR.DF <- t(x)
KonzoData.O.CLR.DF <- KonzoData.O.CLR.DF[4:183,]
rownames(KonzoData.O.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.O.CLR.DF))
KonzoData.O.CLR.DF <- cbind(KonzoData.O.CLR.DF, as.data.frame(KonzoData.O.f@sam_data$Status))
colnames(KonzoData.O.CLR.DF)[colnames(KonzoData.O.CLR.DF)=="KonzoData.O.f@sam_data$Status"] <- "Status"
write.csv(KonzoData.O.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Order_KonzoData_Aldex_MedianCLRValues.csv")  

conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  


#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.O <-  prune_samples((KonzoData.O@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.O@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.O)
KinMasLPZ.O.f <- prune_taxa(f_0.0001, KinMasLPZ.O)   

O.aldex.DF <- as.data.frame(KinMasLPZ.O.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.O <-  prune_samples((KonzoData.O@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.O@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.O)
KinMasHPZ.O.f <- prune_taxa(f_0.0001, KinMasHPZ.O)   

O.aldex.DF <- as.data.frame(KinMasHPZ.O.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_HPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.O <-  prune_samples((KonzoData.O@sam_data$Status == "Kinshasa") | (KonzoData.O@sam_data$Status == "Masimanimba"),  KonzoData.O)
KinMas.O.f <- prune_taxa(f_0.0001, KinMas.O)   

O.aldex.DF <- as.data.frame(KinMas.O.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Kin_vs_Mas_Aldex.csv")  

#Kin vs ULPZ
KinULPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Kinshasa") | (KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.O)
KinULPZ.O.f <- prune_taxa(f_0.0001, KinULPZ.O)   

O.aldex.DF <- as.data.frame(KinULPZ.O.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs ULPZ
MasULPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Masimanimba") | (KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.O)
MasULPZ.O.f <- prune_taxa(f_0.0001, MasULPZ.O)   

O.aldex.DF <- as.data.frame(MasULPZ.O.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs UHPZ
KinUHPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Kinshasa") | (KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.O)
KinUHPZ.O.f <- prune_taxa(f_0.0001, KinUHPZ.O)   

O.aldex.DF <- as.data.frame(KinUHPZ.O.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs UHPZ
MasUHPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Masimanimba") | (KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.O)
MasUHPZ.O.f <- prune_taxa(f_0.0001, MasUHPZ.O)   

O.aldex.DF <- as.data.frame(MasUHPZ.O.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Kinshasa") | (KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.O)
KinKLPZ.O.f <- prune_taxa(f_0.0001, KinKLPZ.O)   

O.aldex.DF <- as.data.frame(KinKLPZ.O.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Kin_vs_KLPZ_Aldex.csv")  

#Mas vs KLPZ
MasKLPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Masimanimba") | (KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.O)
MasKLPZ.O.f <- prune_taxa(f_0.0001, MasKLPZ.O)   

O.aldex.DF <- as.data.frame(MasKLPZ.O.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs KHPZ
KinKHPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Kinshasa") | (KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.O)
KinKHPZ.O.f <- prune_taxa(f_0.0001, KinKHPZ.O)   

O.aldex.DF <- as.data.frame(KinKHPZ.O.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs KHPZ
MasKHPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Masimanimba") | (KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.O)
MasKHPZ.O.f <- prune_taxa(f_0.0001, MasKHPZ.O)   

O.aldex.DF <- as.data.frame(MasKHPZ.O.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_Mas_vs_KHPZ_Aldex.csv")  

#ULPZ vs UHPZ

Control.O <-  prune_samples((KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.O)
Control.O.f <- prune_taxa(f_0.0001, Control.O)   

O.aldex.DF <- as.data.frame(Control.O.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_ULPZ_vs_UHPZ_Aldex.csv")  

#KLPZ vs. KHPZ
Disease.O <-  prune_samples((KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.O)
Disease.O.f <- prune_taxa(f_0.0001, Disease.O)   

O.aldex.DF <- as.data.frame(Disease.O.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30)) #Always Check this
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_KLPZ_vs_KHPZ_Aldex.csv")  

#ULPZ vs. KLPZ

LPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.O)
LPZ.O.f <- prune_taxa(f_0.0001, LPZ.O)   

O.aldex.DF <- as.data.frame(LPZ.O.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_ULPZ_vs_KLPZ_Aldex.csv")  

#UHPZ vs. KHPZ

HPZ.O <-  prune_samples((KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.O)
HPZ.O.f <- prune_taxa(f_0.0001, HPZ.O)   

O.aldex.DF <- as.data.frame(HPZ.O.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30)) #Always Check this
x <- aldex(O.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Order_UHPZ_vs_KHPZ_Aldex.csv")  
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Order_Filtered_Aldex.csv")  
                                     
##FAMILY
x <- read.csv("Kinshasa_Konzo3_Family_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KonzoData
F.aldex.DF <- as.data.frame(KonzoData.F.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.F.Aldex.DF <- x
KonzoData.F.CLR.DF <- t(x)
KonzoData.F.CLR.DF <- KonzoData.F.CLR.DF[4:183,]
rownames(KonzoData.F.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.F.CLR.DF))
KonzoData.F.CLR.DF <- cbind(KonzoData.F.CLR.DF, as.data.frame(KonzoData.F.f@sam_data$Status))
colnames(KonzoData.F.CLR.DF)[colnames(KonzoData.F.CLR.DF)=="KonzoData.F.f@sam_data$Status"] <- "Status"
write.csv(KonzoData.F.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Family_KonzoData_Aldex_MedianCLRValues.csv")  

conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  


#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.F <-  prune_samples((KonzoData.F@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.F@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.F)
KinMasLPZ.F.f <- prune_taxa(f_0.0001, KinMasLPZ.F)   

F.aldex.DF <- as.data.frame(KinMasLPZ.F.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.F <-  prune_samples((KonzoData.F@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.F@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.F)
KinMasHPZ.F.f <- prune_taxa(f_0.0001, KinMasHPZ.F)   

F.aldex.DF <- as.data.frame(KinMasHPZ.F.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_HPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.F <-  prune_samples((KonzoData.F@sam_data$Status == "Kinshasa") | (KonzoData.F@sam_data$Status == "Masimanimba"),  KonzoData.F)
KinMas.F.f <- prune_taxa(f_0.0001, KinMas.F)   

F.aldex.DF <- as.data.frame(KinMas.F.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Kin_vs_Mas_Aldex.csv")  

#Kin vs ULPZ
KinULPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Kinshasa") | (KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.F)
KinULPZ.F.f <- prune_taxa(f_0.0001, KinULPZ.F)   

F.aldex.DF <- as.data.frame(KinULPZ.F.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs ULPZ
MasULPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Masimanimba") | (KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.F)
MasULPZ.F.f <- prune_taxa(f_0.0001, MasULPZ.F)   

F.aldex.DF <- as.data.frame(MasULPZ.F.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs UHPZ
KinUHPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Kinshasa") | (KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.F)
KinUHPZ.F.f <- prune_taxa(f_0.0001, KinUHPZ.F)   

F.aldex.DF <- as.data.frame(KinUHPZ.F.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs UHPZ
MasUHPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Masimanimba") | (KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.F)
MasUHPZ.F.f <- prune_taxa(f_0.0001, MasUHPZ.F)   

F.aldex.DF <- as.data.frame(MasUHPZ.F.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Kinshasa") | (KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.F)
KinKLPZ.F.f <- prune_taxa(f_0.0001, KinKLPZ.F)   

F.aldex.DF <- as.data.frame(KinKLPZ.F.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Kin_vs_KLPZ_Aldex.csv")  

#Mas vs KLPZ
MasKLPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Masimanimba") | (KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.F)
MasKLPZ.F.f <- prune_taxa(f_0.0001, MasKLPZ.F)   

F.aldex.DF <- as.data.frame(MasKLPZ.F.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs KHPZ
KinKHPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Kinshasa") | (KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.F)
KinKHPZ.F.f <- prune_taxa(f_0.0001, KinKHPZ.F)   

F.aldex.DF <- as.data.frame(KinKHPZ.F.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs KHPZ
MasKHPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Masimanimba") | (KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.F)
MasKHPZ.F.f <- prune_taxa(f_0.0001, MasKHPZ.F)   

F.aldex.DF <- as.data.frame(MasKHPZ.F.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_Mas_vs_KHPZ_Aldex.csv")  

#ULPZ vs UHPZ

Control.F <-  prune_samples((KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.F)
Control.F.f <- prune_taxa(f_0.0001, Control.F)   

F.aldex.DF <- as.data.frame(Control.F.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_ULPZ_vs_UHPZ_Aldex.csv")  

#KLPZ vs. KHPZ
Disease.F <-  prune_samples((KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.F)
Disease.F.f <- prune_taxa(f_0.0001, Disease.F)   

F.aldex.DF <- as.data.frame(Disease.F.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30)) #Always Check this
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_KLPZ_vs_KHPZ_Aldex.csv")  

#ULPZ vs. KLPZ

LPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.F)
LPZ.F.f <- prune_taxa(f_0.0001, LPZ.F)   

F.aldex.DF <- as.data.frame(LPZ.F.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_ULPZ_vs_KLPZ_Aldex.csv")  

#UHPZ vs. KHPZ

HPZ.F <-  prune_samples((KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.F)
HPZ.F.f <- prune_taxa(f_0.0001, HPZ.F)   

F.aldex.DF <- as.data.frame(HPZ.F.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30)) #Always Check this
x <- aldex(F.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Family_UHPZ_vs_KHPZ_Aldex.csv")  
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Family_Filtered_Aldex.csv")  


### GENUS
x <- read.csv("Kinshasa_Konzo3_Genus_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KonzoData
G.aldex.DF <- as.data.frame(KonzoData.G.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.G.Aldex.DF <- x
KonzoData.G.CLR.DF <- t(x)
KonzoData.G.CLR.DF <- KonzoData.G.CLR.DF[4:183,]
rownames(KonzoData.G.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.G.CLR.DF))
KonzoData.G.CLR.DF <- cbind(KonzoData.G.CLR.DF, as.data.frame(KonzoData.G.f@sam_data$Status))
colnames(KonzoData.G.CLR.DF)[colnames(KonzoData.G.CLR.DF)=="KonzoData.G.f@sam_data$Status"] <- "Status"
write.csv(KonzoData.G.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Genus_KonzoData_Aldex_MedianCLRValues.csv")  

conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  


#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.G <-  prune_samples((KonzoData.G@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.G)
KinMasLPZ.G.f <- prune_taxa(f_0.0001, KinMasLPZ.G)   

G.aldex.DF <- as.data.frame(KinMasLPZ.G.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
KinMasHPZ.G.f <- prune_taxa(f_0.0001, KinMasHPZ.G)   

G.aldex.DF <- as.data.frame(KinMasHPZ.G.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Masimanimba"),  KonzoData.G)
KinMas.G.f <- prune_taxa(f_0.0001, KinMas.G)   

G.aldex.DF <- as.data.frame(KinMas.G.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")

AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Kin_vs_Mas_Aldex.csv")  



a <- aldex.clr(G.aldex.DF, conds, mc.samples=128, denom="all", verbose=F)
a.tt <- aldex.ttest(a, paired.test=FALSE, verbose=FALSE)
a.effect <- aldex.effect(a, CI=T, verbose=FALSE)
a.all <- data.frame(a.tt,a.effect)
par(mfrow=c(1,2))
aldex.plot(a.all, type="MA", test="wilcoxon")
aldex.plot(x.all, type="MW", test="welch")


#Kin vs ULPZ
KinULPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
KinULPZ.G.f <- prune_taxa(f_0.0001, KinULPZ.G)   

G.aldex.DF <- as.data.frame(KinULPZ.G.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs ULPZ
MasULPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
MasULPZ.G.f <- prune_taxa(f_0.0001, MasULPZ.G)   

G.aldex.DF <- as.data.frame(MasULPZ.G.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs UHPZ
KinUHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
KinUHPZ.G.f <- prune_taxa(f_0.0001, KinUHPZ.G)   

G.aldex.DF <- as.data.frame(KinUHPZ.G.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs UHPZ
MasUHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
MasUHPZ.G.f <- prune_taxa(f_0.0001, MasUHPZ.G)   

G.aldex.DF <- as.data.frame(MasUHPZ.G.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
KinKLPZ.G.f <- prune_taxa(f_0.0001, KinKLPZ.G)   

G.aldex.DF <- as.data.frame(KinKLPZ.G.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Kin_vs_KLPZ_Aldex.csv")  

#Mas vs KLPZ
MasKLPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
MasKLPZ.G.f <- prune_taxa(f_0.0001, MasKLPZ.G)   

G.aldex.DF <- as.data.frame(MasKLPZ.G.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs KHPZ
KinKHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
KinKHPZ.G.f <- prune_taxa(f_0.0001, KinKHPZ.G)   

G.aldex.DF <- as.data.frame(KinKHPZ.G.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs KHPZ
MasKHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
MasKHPZ.G.f <- prune_taxa(f_0.0001, MasKHPZ.G)   

G.aldex.DF <- as.data.frame(MasKHPZ.G.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_Mas_vs_KHPZ_Aldex.csv")  

#ULPZ vs UHPZ

Control.G <-  prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
Control.G.f <- prune_taxa(f_0.0001, Control.G)   

G.aldex.DF <- as.data.frame(Control.G.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_ULPZ_vs_UHPZ_Aldex.csv")  


Control.G.Aldex.DF <- x
Control.G.CLR.DF <- t(x)
Control.G.CLR.DF <- Control.G.CLR.DF[4:63,]
rownames(Control.G.CLR.DF) <- sub('^ULPZ_vs_UHPZ_rab.sample.', '', rownames(Control.G.CLR.DF))
Control.G.CLR.DF <- cbind(Control.G.CLR.DF, as.data.frame(Control.G.f@sam_data$Status))
colnames(Control.G.CLR.DF)[colnames(Control.G.CLR.DF)=="Control.G.f@sam_data$Status"] <- "Status"

write.csv(Control.G.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Genus_ULPZ_vs_UHPZ_Aldex_MedianCLRValues.csv") 



#KLPZ vs. KHPZ
Disease.G <-  prune_samples((KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
Disease.G.f <- prune_taxa(f_0.0001, Disease.G)   

G.aldex.DF <- as.data.frame(Disease.G.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30)) #Always Check this
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_KLPZ_vs_KHPZ_Aldex.csv")  

#ULPZ vs. KLPZ

LPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
LPZ.G.f <- prune_taxa(f_0.0001, LPZ.G)   

G.aldex.DF <- as.data.frame(LPZ.G.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30)) #Always Check this
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_ULPZ_vs_KLPZ_Aldex.csv")  

#UHPZ vs. KHPZ

HPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
HPZ.G.f <- prune_taxa(f_0.0001, HPZ.G)   

G.aldex.DF <- as.data.frame(HPZ.G.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30)) #Always Check this
x <- aldex(G.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Genus_UHPZ_vs_KHPZ_Aldex.csv")  
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Genus_Filtered_Aldex.csv")  




###SPECIES

x <- read.csv("Kinshasa_Konzo3_Species_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KonzoData
S.aldex.DF <- as.data.frame(KonzoData.S.f@otu_table)
conds <- c(rep("Urban", 30), rep("Rural", 150))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Urban_vs_Rural", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_KonzoData_Urban_vs_Rural_Aldex.csv")  

KonzoData.Aldex.DF <- x
KonzoData.CLR.DF <- t(x)
KonzoData.CLR.DF <- KonzoData.CLR.DF[4:183,]
rownames(KonzoData.CLR.DF) <- sub('^Urban_vs_Rural_rab.sample.', '', rownames(KonzoData.CLR.DF))
KonzoData.CLR.DF <- cbind(KonzoData.CLR.DF, as.data.frame(KonzoData.S.f@sam_data$Status))
colnames(KonzoData.CLR.DF)[colnames(KonzoData.CLR.DF)=="KonzoData.S.f@sam_data$Status"] <- "Status"

write.csv(KonzoData.CLR.DF, file = "Kinshasa_Konzo3_Bacteria_Species_KonzoData_Aldex_MedianCLRValues.csv") 

#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))


#sig lactis
#my_comparisons <- list(c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"))

#my_comparisons <- list(c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"))

my_comparisons <- list(c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"),
                         c("Kinshasa", "Konzo_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"))



lact <- ggplot(KonzoData.CLR.DF,aes(x = Status,y = `Lactococcus lactis`)) + 
  geom_boxplot(aes(fill = Status), outlier.size = 0.4, fatten = 0.5) + theme_classic() + ylab("median CLR value") + stat_boxplot(geom ='errorbar') + labs(title = "Lactococcus lactis")
#lact <- lact + geom_jitter(position=position_jitter(0.2), size = 0.2)
lact <- lact + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + theme(plot.title = element_text(face = "italic", size = 8)) + scale_fill_manual(values = konzo_color)
lact <- lact + theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7), axis.title.x = element_blank())
lact <- lact + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)

#sig plant
#my_comparisons <- list(c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"))
                       
plant <- ggplot(KonzoData.CLR.DF,aes(x = Status,y = `Lactobacillus plantarum`)) + 
  geom_boxplot(aes(fill = Status), outlier.size = 0.4, fatten = 0.5) + theme_classic() + ylab("median CLR value") + stat_boxplot(geom ='errorbar') + labs(title = "Lactobacillus plantarum")
plant <- plant + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + theme(plot.title = element_text(face = "italic", size = 8)) + scale_fill_manual(values = konzo_color)
plant <- plant + theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7), axis.title.x = element_blank())
plant <- plant + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)

#sig mesen
#my_comparisons <- list( c("Kinshasa", "Unaffected_High_Prevalence_Zone")) 

mesen <- ggplot(KonzoData.CLR.DF,aes(x = Status,y = `Leuconostoc mesenteroides`)) + 
  geom_boxplot(aes(fill = Status), outlier.size = 0.4, fatten = 0.5) + theme_classic() + ylab("median CLR value") + stat_boxplot(geom ='errorbar') + labs(title = "Leuconostoc mesenteroides")
mesen <- mesen + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + theme(plot.title = element_text(face = "italic", size = 8)) + scale_fill_manual(values = konzo_color)
mesen <- mesen + theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7), axis.title.x = element_blank())
mesen <- mesen + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)


lab_clr <- ggarrange(mesen, lact, plant, labels = c("A","B", "C"), ncol = 3, nrow = 1, font.label = list(size = 7), align = "hv")

tiff(filename = "Kinshasa_Konzo3_LAB_Species_Aldex.tiff", width = 7, height = 3.5, units = "in", res = 600)
lab_clr                                   
dev.off()

#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001



conds <- c(rep("KinMas", 30), rep("Kahemba", 120), rep("KinMas", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_Kahemba", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_KonzoData_KinMas_vs_Kahemba_Aldex.csv")  


#(Kin Mas) vs. (ULPZ KLPZ)
KinMasLPZ.S <-  prune_samples((KonzoData.S@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.S@sam_data$Status != "Konzo_High_Prevalence_Zone"),  KonzoData.S)
KinMasLPZ.S.f <- prune_taxa(f_0.0001, KinMasLPZ.S)   


S.aldex.DF <- as.data.frame(KinMasLPZ.S.f@otu_table)
conds <- c(rep("KinMas", 30), rep("LPZ", 60), rep("KinMas", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_LPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_KonzoData_KinMas_vs_LPZ_Aldex.csv")  


#(Kin Mas) vs. (UHPZ KHPZ)
KinMasHPZ.S <-  prune_samples((KonzoData.S@sam_data$Status != "Unaffected_Low_Prevalence_Zone") & (KonzoData.S@sam_data$Status != "Konzo_Low_Prevalence_Zone"),  KonzoData.S)
KinMasHPZ.S.f <- prune_taxa(f_0.0001, KinMasHPZ.S)   

S.aldex.DF <- as.data.frame(KinMasHPZ.S.f@otu_table)
conds <- c(rep("KinMas", 30), rep("HPZ", 60), rep("KinMas", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KinMas_vs_HPZ", colnames(x), sep = "_")

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_KonzoData_KinMas_vs_HPZ_Aldex.csv")  


#Kin vs Mas
KinMas.S <-  prune_samples((KonzoData.S@sam_data$Status == "Kinshasa") | (KonzoData.S@sam_data$Status == "Masimanimba"),  KonzoData.S)
KinMas.S.f <- prune_taxa(f_0.0001, KinMas.S)   

S.aldex.DF <- as.data.frame(KinMas.S.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")
AD <- x[,70:71]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Kin_vs_Mas_Aldex.csv")  

#KinMas.S.tr <-  prune_samples((KonzoData.S.tr@sam_data$Status == "Kinshasa") | (KonzoData.S.tr@sam_data$Status == "Masimanimba"),  KonzoData.S.tr)
#KinMas.S.tr.f <- prune_taxa(f_0.0001, KinMas.S.tr)   

#S.aldex.DF <- as.data.frame(KinMas.S.f@otu_table)
#conds <- c(rep("Kinshasa", 30), rep("Masimanimba", 30))


#y <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
#           test="t", effect=FALSE)

#colnames(x) <- paste("Kin_vs_Mas", colnames(x), sep = "_")
#AD <- y


#Kin vs ULPZ
KinULPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Kinshasa") | (KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.S)
KinULPZ.S.f <- prune_taxa(f_0.0001, KinULPZ.S)   

S.aldex.DF <- as.data.frame(KinULPZ.S.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("ULPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_ULPZ", colnames(x), sep = "_")

AD <- merge(AD,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD2 <- AD[,-1]
rownames(AD2) <- AD[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Kin_vs_ULPZ_Aldex.csv")  

#Mas vs. ULPZ
MasULPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Masimanimba") | (KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.S)
MasULPZ.S.f <- prune_taxa(f_0.0001, MasULPZ.S)   

S.aldex.DF <- as.data.frame(MasULPZ.S.f@otu_table)
conds <- c(rep("ULPZ", 30), rep("Masimanimba", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_ULPZ", colnames(x), sep = "_")

AD2 <- merge(AD2,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD3 <- AD2[,-1]
rownames(AD3) <- AD2[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Mas_vs_ULPZ_Aldex.csv")  

#Kin vs. UHPZ
KinUHPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Kinshasa") | (KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.S)
KinUHPZ.S.f <- prune_taxa(f_0.0001, KinUHPZ.S)   

S.aldex.DF <- as.data.frame(KinUHPZ.S.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("UHPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_UHPZ", colnames(x), sep = "_")

AD3 <- merge(AD3,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD4 <- AD3[,-1]
rownames(AD4) <- AD3[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Kin_vs_UHPZ_Aldex.csv")  

#Mas vs. UHPZ
MasUHPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Masimanimba") | (KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.S)
MasUHPZ.S.f <- prune_taxa(f_0.0001, MasUHPZ.S)   

S.aldex.DF <- as.data.frame(MasUHPZ.S.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("Masimanimba", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_UHPZ", colnames(x), sep = "_")

AD4 <- merge(AD4,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD5 <- AD4[,-1]
rownames(AD5) <- AD4[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Mas_vs_UHPZ_Aldex.csv")  

#Kin vs KLPZ
KinKLPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Kinshasa") | (KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.S)
KinKLPZ.S.f <- prune_taxa(f_0.0001, KinKLPZ.S)   

S.aldex.DF <- as.data.frame(KinKLPZ.S.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KLPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",  effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KLPZ", colnames(x), sep = "_")

AD5 <- merge(AD5,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD6 <- AD5[,-1]
rownames(AD6) <- AD5[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Kin_vs_KLPZ_Aldex.csv")  
 
                                     
#Mas vs. KLPZ
MasKLPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Masimanimba") | (KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.S)
MasKLPZ.S.f <- prune_taxa(f_0.0001, MasKLPZ.S)   

S.aldex.DF <- as.data.frame(MasKLPZ.S.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("Masimanimba", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KLPZ", colnames(x), sep = "_")

AD6 <- merge(AD6,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD7 <- AD6[,-1]
rownames(AD7) <- AD6[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Mas_vs_KLPZ_Aldex.csv")  

#Kin vs. KHPZ
KinKHPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Kinshasa") | (KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.S)
KinKHPZ.S.f <- prune_taxa(f_0.0001, KinKHPZ.S)   

S.aldex.DF <- as.data.frame(KinKHPZ.S.f@otu_table)
conds <- c(rep("Kinshasa", 30), rep("KHPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Kin_vs_KHPZ", colnames(x), sep = "_")

AD7 <- merge(AD7,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD8 <- AD7[,-1]
rownames(AD8) <- AD7[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Kin_vs_KHPZ_Aldex.csv")  

#Mas vs. KHPZ
MasKHPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Masimanimba") | (KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.S)
MasKHPZ.S.f <- prune_taxa(f_0.0001, MasKHPZ.S)   

S.aldex.DF <- as.data.frame(MasKHPZ.S.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("Masimanimba", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("Mas_vs_KHPZ", colnames(x), sep = "_")

AD8 <- merge(AD8,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD9 <- AD8[,-1]
rownames(AD9) <- AD8[,1]


write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_Mas_vs_KHPZ_Aldex.csv") 

#ULPZ vs. UHPZ
Control.S <-  prune_samples((KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.S)
Control.S.f <- prune_taxa(f_0.0001, Control.S)   

S.aldex.DF <- as.data.frame(Control.S.f@otu_table)
conds <- c(rep("UHPZ", 30), rep("ULPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_UHPZ", colnames(x), sep = "_")

AD9 <- merge(AD9,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD10 <- AD9[,-1]
rownames(AD10) <- AD9[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_ULPZ_vs_UHPZ_Aldex.csv") 

#KLPZ vs. KHPZ
Disease.S <-  prune_samples((KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.S)
Disease.S.f <- prune_taxa(f_0.0001, Disease.S)   

S.aldex.DF <- as.data.frame(Disease.S.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("KLPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("KLPZ_vs_KHPZ", colnames(x), sep = "_")

AD10 <- merge(AD10,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD11 <- AD10[,-1]
rownames(AD11) <- AD10[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_KLPZ_vs_KHPZ_Aldex.csv") 

#ULPZ vs. KLPZ
LPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.S)
LPZ.S.f <- prune_taxa(f_0.0001, LPZ.S)   

S.aldex.DF <- as.data.frame(LPZ.S.f@otu_table)
conds <- c(rep("KLPZ", 30), rep("ULPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t", effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("ULPZ_vs_KLPZ", colnames(x), sep = "_")

AD11 <- merge(AD11,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD12 <- AD11[,-1]
rownames(AD12) <- AD11[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_ULPZ_vs_KLPZ_Aldex.csv") 

#UHPZ vs. KHPZ
HPZ.S <-  prune_samples((KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.S)
HPZ.S.f <- prune_taxa(f_0.0001, HPZ.S)   

S.aldex.DF <- as.data.frame(HPZ.S.f@otu_table)
conds <- c(rep("KHPZ", 30), rep("UHPZ", 30))
x <- aldex(S.aldex.DF, conds, mc.samples=128, denom="all",
           test="t",effect=TRUE, include.sample.summary=TRUE)

colnames(x) <- paste("UHPZ_vs_KHPZ", colnames(x), sep = "_")

AD12 <- merge(AD12,x[,70:71],by='row.names',all=TRUE, sort = FALSE)
AD13 <- AD12[,-1]
rownames(AD13) <- AD12[,1]

write.csv(x, file = "Kinshasa_Konzo3_Bacteria_Species_UHPZ_vs_KHPZ_Aldex.csv") 
write.csv(AD13, file = "Kinshasa_Konzo3_Bacteria_Species_Filtered_Aldex.csv") 

                                     
                                     
### MANN WHITNEY_WILCOX TEST (with BH correction)
#Supplemental File 4 where the saved WT (results from the mann whitney test) are joined into one excel sheet for all the different comparisions, and each tab is each taxa rank                           
                           
#Bacteria Phylum
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Phylum")
x <- read.csv("Kinshasa_Konzo3_Phylum_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                                                                              
#KINSHASA AND MASIMANIMBA
KinMas.P <-  prune_samples(KonzoData.P@sam_data$Status == "Kinshasa" | KonzoData.P@sam_data$Status == "Masimanimba", KonzoData.P)
KinMas.P.tr <-  transform_sample_counts(KinMas.P, function(x) x / sum(x))
KinMas.P.tr.f <-  prune_taxa(f_0.0001, KinMas.P.tr)

#MWW 
                                               
P <- KinMas.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- KinMas.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                          
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Kinshasa vs. Masi-manimba p-value" ,"Kinshasa vs. Masi-manimba p-value adjusted")
for (i in 1:(ncol(P.tr.DF)-1)) #ncol - 1 because the last column is Status, so need to run test on that
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinMas_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.P.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.P.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.P.tr.f)                                        
                                        
write.csv(KinMas.P.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.P.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.P.tr.f.status <- merge_samples(KinMas.P.tr.f, KinMas.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.P.tr.f.status <- transform_sample_counts(KinMas.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.P.tr.f.status)                                        
KinMas.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.P.tr.f.status)                                        
                                                                                                
write.csv(t(KinMas.P.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.P.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                                             
MWW_phylum <- WT
                                                                                                                         
#KINSHASA AND UNAFFECTED LPZ
                                             
KinULPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Kinshasa" | KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.P)
KinULPZ.P.tr <-  transform_sample_counts(KinULPZ.P, function(x) x / sum(x))
KinULPZ.P.tr.f <-  prune_taxa(f_0.0001, KinULPZ.P.tr)
#MWW 
                                               
P <- KinULPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- KinULPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Kinshasa vs. ULPZ p-value",  "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "KinULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#KinULPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinULPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.P.tr.f)                                        
                                        
#write.csv(KinULPZ.P.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinULPZ.P.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.P.tr.f.status <- merge_samples(KinULPZ.P.tr.f, KinULPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.P.tr.f.status <- transform_sample_counts(KinULPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#KinULPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.P.tr.f.status)                                        
#KinULPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.P.tr.f.status)                                        
                                                                                                
#write.csv(t(KinULPZ.P.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinULPZ.P.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)
 
#MASIMANIMBA AND UNAFFECTED LPZ                                        
MasULPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba" | KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.P)
MasULPZ.P.tr <- transform_sample_counts(MasULPZ.P, function(x) x / sum(x)) 
MasULPZ.P.tr.f <- prune_taxa(f_0.0001, MasULPZ.P.tr)
                                       
P <- MasULPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- MasULPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Masi-manimba vs. ULPZ p-value", "Masi-manimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}

WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.P.tr.f)                                        

write.csv(MasULPZ.P.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.P.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

MasULPZ.P.tr.f.status <- merge_samples(MasULPZ.P.tr.f, MasULPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.P.tr.f.status <- transform_sample_counts(MasULPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.P.tr.f.status)                                        
MasULPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.P.tr.f.status)                                        

write.csv(t(MasULPZ.P.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.P.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)

                                        
#KINSHASA AND UNAFFECTED HPZ                                            
KinUHPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Kinshasa" | KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.P)
KinUHPZ.P.tr <-  transform_sample_counts(KinUHPZ.P, function(x) x / sum(x))
KinUHPZ.P.tr.f <-  prune_taxa(f_0.0001, KinUHPZ.P.tr)
#MWW 
                                               
P <- KinUHPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- KinUHPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Kinshasa vs. UHPZ p-value",  "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.P.tr.f)                                        
                                        
write.csv(KinUHPZ.P.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.P.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.P.tr.f.status <- merge_samples(KinUHPZ.P.tr.f, KinUHPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.P.tr.f.status <- transform_sample_counts(KinUHPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.P.tr.f.status)                                        
KinUHPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.P.tr.f.status)                                        
                                                                                                
write.csv(t(KinUHPZ.P.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.P.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)                                             
                                                                                               
#MASIMANIMBA AND UNAFFECTED HPZ  
MasUHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba" | KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.P)
MasUHPZ.P.tr <- transform_sample_counts(MasUHPZ.P, function(x) x / sum(x)) 
MasUHPZ.P.tr.f <- prune_taxa(f_0.0001, MasUHPZ.P.tr)
                                       
P <- MasUHPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- MasUHPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Masi-manimba vs. UHPZ p-value", "Masi-manimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}

WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "MasUHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasUHPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#MasUHPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.P.tr.f)                                        

write.csv(MasUHPZ.P.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(MasUHPZ.P.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

MasUHPZ.P.tr.f.status <- merge_samples(MasUHPZ.P.tr.f, MasUHPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.P.tr.f.status <- transform_sample_counts(MasUHPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasUHPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.P.tr.f.status)                                        
#MasUHPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.P.tr.f.status)                                        

write.csv(t(MasUHPZ.P.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(MasUHPZ.P.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)                                                 
                                                 
#KINSHASA AND KONZO LPZ
                                             
KinKLPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Kinshasa" | KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.P)
KinKLPZ.P.tr <-  transform_sample_counts(KinKLPZ.P, function(x) x / sum(x))
KinKLPZ.P.tr.f <-  prune_taxa(f_0.0001, KinKLPZ.P.tr)
                                         
#MWW 
                                               
P <- KinKLPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- KinKLPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Kinshasa vs. KLPZ p-value",  "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "KinKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#KinKLPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinKLPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.P.tr.f)                                        
                                        
#write.csv(KinKLPZ.P.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinKLPZ.P.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.P.tr.f.status <- merge_samples(KinKLPZ.P.tr.f, KinKLPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.P.tr.f.status <- transform_sample_counts(KinKLPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#KinKLPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.P.tr.f.status)                                        
#KinKLPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.P.tr.f.status)                                        
                                                                                                
#write.csv(t(KinKLPZ.P.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinKLPZ.P.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)                                                 
                                                 
#MASIMANIMBA AND KONZO LPZ
MasKLPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba" | KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.P)
MasKLPZ.P.tr <- transform_sample_counts(MasKLPZ.P, function(x) x / sum(x)) 
MasKLPZ.P.tr.f <- prune_taxa(f_0.0001, MasKLPZ.P.tr)
                                       
P <- MasKLPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- MasKLPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Masi-manimba vs. KLPZ p-value", "Masi-manimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}

WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.P.tr.f)                                        

write.csv(MasKLPZ.P.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.P.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

MasKLPZ.P.tr.f.status <- merge_samples(MasKLPZ.P.tr.f, MasKLPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.P.tr.f.status <- transform_sample_counts(MasKLPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.P.tr.f.status)                                        
MasKLPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.P.tr.f.status)                                        

write.csv(t(MasKLPZ.P.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.P.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE) 
                                                 
#KINSHASA AND KONZO HPZ
KinKHPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Kinshasa" | KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
KinKHPZ.P.tr <-  transform_sample_counts(KinKHPZ.P, function(x) x / sum(x))
KinKHPZ.P.tr.f <-  prune_taxa(f_0.0001, KinKHPZ.P.tr)
                                         
#MWW 
                                               
P <- KinKHPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- KinKHPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Kinshasa vs. KHPZ p-value",  "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.P.tr.f)                                        
                                        
write.csv(KinKHPZ.P.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.P.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.P.tr.f.status <- merge_samples(KinKHPZ.P.tr.f, KinKHPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.P.tr.f.status <- transform_sample_counts(KinKHPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.P.tr.f.status)                                        
KinKHPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.P.tr.f.status)                                        
                                                                                                
write.csv(t(KinKHPZ.P.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.P.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)                                                  
                                                                                                
#MASIMANIMBA AND KONZO HPZ
MasKHPZ.P <- prune_samples(KonzoData.P@sam_data$Status == "Masimanimba" | KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
MasKHPZ.P.tr <- transform_sample_counts(MasKHPZ.P, function(x) x / sum(x)) 
MasKHPZ.P.tr.f <- prune_taxa(f_0.0001, MasKHPZ.P.tr)
                                       
P <- MasKHPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- MasKHPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "Masi-manimba vs. KHPZ p-value", "Masi-manimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}

WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.P.tr.f)                                        

write.csv(MasKHPZ.P.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.P.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

MasKHPZ.P.tr.f.status <- merge_samples(MasKHPZ.P.tr.f, MasKHPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.P.tr.f.status <- transform_sample_counts(MasKHPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.P.tr.f.status)                                        
MasKHPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.P.tr.f.status)                                        

write.csv(t(MasKHPZ.P.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.P.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)   
                                                 
#CONTROL (UNAFFECTED)
Control.P <-  prune_samples(KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.P)
Control.P.tr <-  transform_sample_counts(Control.P, function(x) x / sum(x))
Control.P.tr.f <-  prune_taxa(f_0.0001, Control.P.tr)
                                         
#MWW 
                                               
P <- Control.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- Control.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "ULPZ vs. UHPZ p-value",  "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Control_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.P.tr.f.0.05 <- prune_taxa(ls_0.05,Control.P.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.P.tr.f.0.01 <- prune_taxa(ls_0.01,Control.P.tr.f)                                        
                                        
write.csv(Control.P.tr.f.0.05@otu_table, file = "./Control_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.P.tr.f.0.01@otu_table, file = "./Control_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Control.P.tr.f.status <- merge_samples(Control.P.tr.f, Control.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.P.tr.f.status <- transform_sample_counts(Control.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.P.tr.f.status)                                        
Control.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.P.tr.f.status)                                        
                                                                                                
write.csv(t(Control.P.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.P.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)
                                                 
#DISEASE (KONZO)
Disease.P <-  prune_samples(KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
Disease.P.tr <-  transform_sample_counts(Disease.P, function(x) x / sum(x))
Disease.P.tr.f <-  prune_taxa(f_0.0001, Disease.P.tr)
                                         
#MWW 
                                               
P <- Disease.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- Disease.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "KLPZ vs. KHPZ p-value",  "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Disease_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "Disease_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "Disease_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#Disease.P.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#Disease.P.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.P.tr.f)                                        
                                        
#write.csv(Disease.P.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.P.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Disease.P.tr.f.status <- merge_samples(Disease.P.tr.f, Disease.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Disease.P.tr.f.status <- transform_sample_counts(Disease.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#Disease.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.P.tr.f.status)                                        
#Disease.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.P.tr.f.status)                                        
                                                                                                
#write.csv(t(Disease.P.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.P.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)

#LPZ
LPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.P@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.P)
LPZ.P.tr <-  transform_sample_counts(LPZ.P, function(x) x / sum(x))
LPZ.P.tr.f <-  prune_taxa(f_0.0001, LPZ.P.tr)
                                         
#MWW 
                                               
P <- LPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- LPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "ULPZ vs. KLPZ p-value",  "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "LPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "LPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "LPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#LPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#LPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.P.tr.f)                                        
                                        
#write.csv(LPZ.P.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.P.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
LPZ.P.tr.f.status <- merge_samples(LPZ.P.tr.f, LPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
LPZ.P.tr.f.status <- transform_sample_counts(LPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.P.tr.f.status)                                        
#LPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.P.tr.f.status)                                        
                                                                                                
#write.csv(t(LPZ.P.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.P.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)  
                                             
#HPZ
HPZ.P <-  prune_samples(KonzoData.P@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData.P@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.P)
HPZ.P.tr <-  transform_sample_counts(HPZ.P, function(x) x / sum(x))
HPZ.P.tr.f <-  prune_taxa(f_0.0001, HPZ.P.tr)
                                         
#MWW 
                                               
P <- HPZ.P.tr.f
                                               
P.tr_META <- as.data.frame(P@sam_data)
P.tr_OTU <- as.data.frame(t(P@otu_table))
P.tr.DF <- cbind(P.tr_OTU, P.tr_META$Status)

colnames(P.tr.DF)[colnames(P.tr.DF)=="P.tr_META$Status"] <- "Status"
for (i in 1:nrow(P.tr.DF))
  {P.tr.DF[i,]$Status <- HPZ.P.tr.f@sam_data[rownames(P.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(P.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Phylum", "UHPZ vs. KHPZ p-value",  "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(P.tr.DF)-1))
{
  wt <- wilcox.test(P.tr.DF[,i] ~P.tr.DF$Status, data = P.tr.DF)
  WT[i,1] = colnames(P.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "HPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "HPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "HPZ_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#HPZ.P.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.P.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#HPZ.P.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.P.tr.f)                                        
                                        
#write.csv(HPZ.P.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.P.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Phylum_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
HPZ.P.tr.f.status <- merge_samples(HPZ.P.tr.f, HPZ.P.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
HPZ.P.tr.f.status <- transform_sample_counts(HPZ.P.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.P.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.P.tr.f.status)                                        
#HPZ.P.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.P.tr.f.status)                                        
                                                                                                
#write.csv(t(HPZ.P.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.P.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Phylum_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                             
MWW_phylum <- merge(MWW_phylum,WT,by="Bacteria Phylum", sort = FALSE)  
write.csv(MWW_phylum, file = "Kinshasa_Konzo3_Bacteria_Phylum_f_0.0001_ByStatus_WilcoxTest_BH.csv")                                            
                                             
#Bacteria Class
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Class")
x <- read.csv("Kinshasa_Konzo3_Class_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                             
#KINSHASA AND MASIMANIMBA
KinMas.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa" | KonzoData.C@sam_data$Status == "Masimanimba", KonzoData.C)                                        
KinMas.C.tr <- transform_sample_counts(KinMas.C, function(x) x / sum(x))                                             
KinMas.C.tr.f <- prune_taxa(f_0.0001, KinMas.C.tr)  

C <- KinMas.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- KinMas.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Kinshasa vs. Masi-manimba p-value", "Kinshasa vs. Masi-manimba p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinMas_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.C.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.C.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.C.tr.f)                                        
                                        
write.csv(KinMas.C.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.C.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.C.tr.f.status <- merge_samples(KinMas.C.tr.f, KinMas.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.C.tr.f.status <- transform_sample_counts(KinMas.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.C.tr.f.status)                                        
KinMas.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.C.tr.f.status)                                        
                                                                                                
write.csv(t(KinMas.C.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.C.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                                                       
MWW_class <- WT
                                        
#KINSHASA AND UNAFFECTED LPZ
KinULPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa" | KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.C)                                        
KinULPZ.C.tr <- transform_sample_counts(KinULPZ.C, function(x) x / sum(x))                                             
KinULPZ.C.tr.f <- prune_taxa(f_0.0001, KinULPZ.C.tr)  

C <- KinULPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- KinULPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Kinshasa vs. ULPZ p-value", "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinULPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinULPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.C.tr.f)                                        
                                        
write.csv(KinULPZ.C.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinULPZ.C.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.C.tr.f.status <- merge_samples(KinULPZ.C.tr.f, KinULPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.C.tr.f.status <- transform_sample_counts(KinULPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinULPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.C.tr.f.status)                                        
KinULPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(KinULPZ.C.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinULPZ.C.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)
                                                
                                                
#MASIMANIMBA AND UNAFFECTED LPZ
MasULPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba" | KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.C)                                        
MasULPZ.C.tr <- transform_sample_counts(MasULPZ.C, function(x) x / sum(x))                                             
MasULPZ.C.tr.f <- prune_taxa(f_0.0001, MasULPZ.C.tr)  

C <- MasULPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- MasULPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Masimanimba vs. ULPZ p-value", "Masimanimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.C.tr.f)                                        
                                        
write.csv(MasULPZ.C.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.C.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasULPZ.C.tr.f.status <- merge_samples(MasULPZ.C.tr.f, MasULPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.C.tr.f.status <- transform_sample_counts(MasULPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.C.tr.f.status)                                        
MasULPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(MasULPZ.C.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.C.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)
                                                 
#KINSHASA AND UNAFFECTED HPZ 
KinUHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa" | KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.C)                                        
KinUHPZ.C.tr <- transform_sample_counts(KinUHPZ.C, function(x) x / sum(x))                                             
KinUHPZ.C.tr.f <- prune_taxa(f_0.0001, KinUHPZ.C.tr)  

C <- KinUHPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- KinUHPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Kinshasa vs. UHPZ p-value", "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.C.tr.f)                                        
                                        
write.csv(KinUHPZ.C.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.C.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.C.tr.f.status <- merge_samples(KinUHPZ.C.tr.f, KinUHPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.C.tr.f.status <- transform_sample_counts(KinUHPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.C.tr.f.status)                                        
KinUHPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(KinUHPZ.C.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.C.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)
                                
#MASIMANIMBA AND UNAFFECTED HPZ
MasUHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba" | KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.C)                                        
MasUHPZ.C.tr <- transform_sample_counts(MasUHPZ.C, function(x) x / sum(x))                                             
MasUHPZ.C.tr.f <- prune_taxa(f_0.0001, MasUHPZ.C.tr)  

C <- MasUHPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- MasUHPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Masimanimba vs. UHPZ p-value", "Masimanimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasUHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasUHPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.C.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#MasUHPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.C.tr.f)                                        
                                        
write.csv(MasUHPZ.C.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(MasUHPZ.C.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasUHPZ.C.tr.f.status <- merge_samples(MasUHPZ.C.tr.f, MasUHPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.C.tr.f.status <- transform_sample_counts(MasUHPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasUHPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.C.tr.f.status)                                        
#MasUHPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(MasUHPZ.C.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(MasUHPZ.C.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)
                                                
#KINSHASA AND KONZO LPZ
                                                 
KinKLPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa" | KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.C)                                        
KinKLPZ.C.tr <- transform_sample_counts(KinKLPZ.C, function(x) x / sum(x))                                             
KinKLPZ.C.tr.f <- prune_taxa(f_0.0001, KinKLPZ.C.tr)  

C <- KinKLPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- KinKLPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Kinshasa vs. KLPZ p-value", "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "KinKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#KinKLPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.C.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinKLPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.C.tr.f)                                        
                                        
#write.csv(KinKLPZ.C.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinKLPZ.C.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.C.tr.f.status <- merge_samples(KinKLPZ.C.tr.f, KinKLPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.C.tr.f.status <- transform_sample_counts(KinKLPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#KinKLPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.C.tr.f.status)                                        
#KinKLPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.C.tr.f.status)                                        
                                                                                                
#write.csv(t(KinKLPZ.C.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinKLPZ.C.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)
                                                
                                                                                     
#MASIMANIMBA AND KONZO LPZ
MasKLPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba" | KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.C)                                        
MasKLPZ.C.tr <- transform_sample_counts(MasKLPZ.C, function(x) x / sum(x))                                             
MasKLPZ.C.tr.f <- prune_taxa(f_0.0001, MasKLPZ.C.tr)  

C <- MasKLPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- MasKLPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Masimanimba vs. KLPZ p-value", "Masimanimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.C.tr.f)                                        
                                        
write.csv(MasKLPZ.C.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.C.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKLPZ.C.tr.f.status <- merge_samples(MasKLPZ.C.tr.f, MasKLPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.C.tr.f.status <- transform_sample_counts(MasKLPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.C.tr.f.status)                                        
MasKLPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(MasKLPZ.C.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.C.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)                                                 

#KINSHASA AND KONZO HPZ 
KinKHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Kinshasa" | KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)                                        
KinKHPZ.C.tr <- transform_sample_counts(KinKHPZ.C, function(x) x / sum(x))                                             
KinKHPZ.C.tr.f <- prune_taxa(f_0.0001, KinKHPZ.C.tr)  

C <- KinKHPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- KinKHPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Kinshasa vs. KHPZ p-value", "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.C.tr.f)                                        
                                        
write.csv(KinKHPZ.C.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.C.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.C.tr.f.status <- merge_samples(KinKHPZ.C.tr.f, KinKHPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.C.tr.f.status <- transform_sample_counts(KinKHPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.C.tr.f.status)                                        
KinKHPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(KinKHPZ.C.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.C.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)                                                 
                                
#MASIMANIMBA AND KONZO HPZ                                                
MasKHPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Masimanimba" | KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)                                        
MasKHPZ.C.tr <- transform_sample_counts(MasKHPZ.C, function(x) x / sum(x))                                             
MasKHPZ.C.tr.f <- prune_taxa(f_0.0001, MasKHPZ.C.tr)  

C <- MasKHPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- MasKHPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "Masimanimba vs. KHPZ p-value", "Masimanimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.C.tr.f)                                        
                                        
write.csv(MasKHPZ.C.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.C.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKHPZ.C.tr.f.status <- merge_samples(MasKHPZ.C.tr.f, MasKHPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.C.tr.f.status <- transform_sample_counts(MasKHPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.C.tr.f.status)                                        
MasKHPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.C.tr.f.status)                                        
                                                                                                
write.csv(t(MasKHPZ.C.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.C.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)  
                                                 
#CONTROL (UNAFFECTED LPZ vs. HPZ)
Control.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.C)                                        
Control.C.tr <- transform_sample_counts(Control.C, function(x) x / sum(x))                                             
Control.C.tr.f <- prune_taxa(f_0.0001, Control.C.tr)  

C <- Control.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- Control.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "ULPZ vs. UHPZ p-value", "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Control_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.C.tr.f.0.05 <- prune_taxa(ls_0.05,Control.C.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.C.tr.f.0.01 <- prune_taxa(ls_0.01,Control.C.tr.f)                                        
                                        
write.csv(Control.C.tr.f.0.05@otu_table, file = "./Control_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.C.tr.f.0.01@otu_table, file = "./Control_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Control.C.tr.f.status <- merge_samples(Control.C.tr.f, Control.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.C.tr.f.status <- transform_sample_counts(Control.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.C.tr.f.status)                                        
Control.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.C.tr.f.status)                                        
                                                                                                
write.csv(t(Control.C.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.C.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)   
                                                 
#DISEASE (KONZO LPZ vs. HPZ)
Disease.C <- prune_samples(KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)                                        
Disease.C.tr <- transform_sample_counts(Disease.C, function(x) x / sum(x))                                             
Disease.C.tr.f <- prune_taxa(f_0.0001, Disease.C.tr)  

C <- Disease.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- Disease.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "KLPZ vs. KHPZ p-value", "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Disease_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Disease_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "Disease_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Disease.C.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.C.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#Disease.C.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.C.tr.f)                                        
                                        
write.csv(Disease.C.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.C.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Disease.C.tr.f.status <- merge_samples(Disease.C.tr.f, Disease.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Disease.C.tr.f.status <- transform_sample_counts(Disease.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Disease.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.C.tr.f.status)                                        
#Disease.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.C.tr.f.status)                                        
                                                                                                
write.csv(t(Disease.C.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.C.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)                                                  
                                          
#LPZ (UNAFFECTED LPZ vs. KONZO LPZ)
LPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.C@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.C)                                        
LPZ.C.tr <- transform_sample_counts(LPZ.C, function(x) x / sum(x))                                             
LPZ.C.tr.f <- prune_taxa(f_0.0001, LPZ.C.tr)  

C <- LPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- LPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "ULPZ vs. KLPZ p-value", "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "LPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
#WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "LPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "LPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#LPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.C.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#LPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.C.tr.f)                                        
                                        
#write.csv(LPZ.C.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.C.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
LPZ.C.tr.f.status <- merge_samples(LPZ.C.tr.f, LPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
LPZ.C.tr.f.status <- transform_sample_counts(LPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.C.tr.f.status)                                        
#LPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.C.tr.f.status)                                        
                                                                                                
#write.csv(t(LPZ.C.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.C.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                  
                                                
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)                                                   
                                      
#HPZ (UNAFFECTED HPZ vs. KONZO HPZ)
HPZ.C <- prune_samples(KonzoData.C@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData.C@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.C)                                        
HPZ.C.tr <- transform_sample_counts(HPZ.C, function(x) x / sum(x))                                             
HPZ.C.tr.f <- prune_taxa(f_0.0001, HPZ.C.tr)  

C <- HPZ.C.tr.f
                                               
C.tr_META <- as.data.frame(C@sam_data)
C.tr_OTU <- as.data.frame(t(C@otu_table))
C.tr.DF <- cbind(C.tr_OTU, C.tr_META$Status)

colnames(C.tr.DF)[colnames(C.tr.DF)=="C.tr_META$Status"] <- "Status"
for (i in 1:nrow(C.tr.DF))
  {C.tr.DF[i,]$Status <- HPZ.C.tr.f@sam_data[rownames(C.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(C.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Class", "UHPZ vs. KHPZ p-value", "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(C.tr.DF)-1))
{
  wt <- wilcox.test(C.tr.DF[,i] ~C.tr.DF$Status, data = C.tr.DF)
  WT[i,1] = colnames(C.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "HPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
#WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "HPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "HPZ_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#HPZ.C.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.C.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#HPZ.C.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.C.tr.f)                                        
                                        
#write.csv(HPZ.C.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.C.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Class_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
HPZ.C.tr.f.status <- merge_samples(HPZ.C.tr.f, HPZ.C.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
HPZ.C.tr.f.status <- transform_sample_counts(HPZ.C.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.C.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.C.tr.f.status)                                        
#HPZ.C.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.C.tr.f.status)                                        
                                                                                                
#write.csv(t(HPZ.C.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.C.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Class_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_class <- merge(MWW_class,WT,by="Bacteria Class", sort = FALSE)  
write.csv(MWW_class, file = "Kinshasa_Konzo3_Bacteria_Class_f_0.0001_ByStatus_WilcoxTest_BH.csv")                                            
                                                 
#Bacteria Order
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Order")
x <- read.csv("Kinshasa_Konzo3_Order_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#KINSHASA AND MASIMANIMBA
KinMas.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa" | KonzoData.O@sam_data$Status == "Masimanimba", KonzoData.O)                                        
KinMas.O.tr <- transform_sample_counts(KinMas.O, function(x) x / sum(x))                                             
KinMas.O.tr.f <- prune_taxa(f_0.0001, KinMas.O.tr)  

O <- KinMas.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- KinMas.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Kinshasa vs. Masimanimba p-value", "Kinshasa vs. Masimanimba p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinMas_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.O.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.O.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.O.tr.f)                                        
                                        
write.csv(KinMas.O.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.O.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.O.tr.f.status <- merge_samples(KinMas.O.tr.f, KinMas.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.O.tr.f.status <- transform_sample_counts(KinMas.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.O.tr.f.status)                                        
KinMas.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.O.tr.f.status)                                        
                                                                                                
write.csv(t(KinMas.O.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.O.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_order <- WT
                                             
#KINSHASA AND UNAFFECTED LPZ
KinULPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa" | KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.O)                                        
KinULPZ.O.tr <- transform_sample_counts(KinULPZ.O, function(x) x / sum(x))                                             
KinULPZ.O.tr.f <- prune_taxa(f_0.0001, KinULPZ.O.tr)  

O <- KinULPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- KinULPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Kinshasa vs. ULPZ p-value", "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinULPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinULPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.O.tr.f)                                        
                                        
write.csv(KinULPZ.O.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinULPZ.O.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.O.tr.f.status <- merge_samples(KinULPZ.O.tr.f, KinULPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.O.tr.f.status <- transform_sample_counts(KinULPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinULPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.O.tr.f.status)                                        
#KinULPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(KinULPZ.O.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinULPZ.O.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                
#MASIMANIMBA AND UNAFFECTED LPZ  
MasULPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba" | KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.O)                                        
MasULPZ.O.tr <- transform_sample_counts(MasULPZ.O, function(x) x / sum(x))                                             
MasULPZ.O.tr.f <- prune_taxa(f_0.0001, MasULPZ.O.tr)  

O <- MasULPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- MasULPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Masimanimba vs. ULPZ p-value", "Masimanimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.O.tr.f)                                        
                                        
write.csv(MasULPZ.O.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.O.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasULPZ.O.tr.f.status <- merge_samples(MasULPZ.O.tr.f, MasULPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.O.tr.f.status <- transform_sample_counts(MasULPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.O.tr.f.status)                                        
MasULPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(MasULPZ.O.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.O.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                 

#KINSHASA vs UNAFFECTED HPZ
KinUHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa" | KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.O)                                        
KinUHPZ.O.tr <- transform_sample_counts(KinUHPZ.O, function(x) x / sum(x))                                             
KinUHPZ.O.tr.f <- prune_taxa(f_0.0001, KinUHPZ.O.tr)  

O <- KinUHPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- KinUHPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Kinshasa vs. UHPZ p-value", "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.O.tr.f)                                        
                                        
write.csv(KinUHPZ.O.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.O.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.O.tr.f.status <- merge_samples(KinUHPZ.O.tr.f, KinUHPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.O.tr.f.status <- transform_sample_counts(KinUHPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.O.tr.f.status)                                        
KinUHPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(KinUHPZ.O.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.O.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                
#MASIMANIMBA vs. UNAFFECTED HPZ  
MasUHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba" | KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.O)                                        
MasUHPZ.O.tr <- transform_sample_counts(MasUHPZ.O, function(x) x / sum(x))                                             
MasUHPZ.O.tr.f <- prune_taxa(f_0.0001, MasUHPZ.O.tr)  

O <- MasUHPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- MasUHPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Masimanimba vs. UHPZ p-value", "Masimanimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "MasUHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasUHPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#MasUHPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.O.tr.f)                                        
                                        
write.csv(MasUHPZ.O.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(MasUHPZ.O.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasUHPZ.O.tr.f.status <- merge_samples(MasUHPZ.O.tr.f, MasUHPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.O.tr.f.status <- transform_sample_counts(MasUHPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasUHPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.O.tr.f.status)                                        
#MasUHPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(MasUHPZ.O.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(MasUHPZ.O.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                 
#KINSHASA AND KONZO LPZ
KinKLPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa" | KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.O)                                        
KinKLPZ.O.tr <- transform_sample_counts(KinKLPZ.O, function(x) x / sum(x))                                             
KinKLPZ.O.tr.f <- prune_taxa(f_0.0001, KinKLPZ.O.tr)  

O <- KinKLPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- KinKLPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Kinshasa vs. KLPZ p-value", "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKLPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinKLPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.O.tr.f)                                        
                                        
write.csv(KinKLPZ.O.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinKLPZ.O.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.O.tr.f.status <- merge_samples(KinKLPZ.O.tr.f, KinKLPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.O.tr.f.status <- transform_sample_counts(KinKLPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKLPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.O.tr.f.status)                                        
#KinKLPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(KinKLPZ.O.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinKLPZ.O.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE) 
                                                 
#MASIMANIMBA AND KONZO LPZ                                        
MasKLPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba" | KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.O)                                        
MasKLPZ.O.tr <- transform_sample_counts(MasKLPZ.O, function(x) x / sum(x))                                             
MasKLPZ.O.tr.f <- prune_taxa(f_0.0001, MasKLPZ.O.tr)  

O <- MasKLPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- MasKLPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Masimanimba vs. KLPZ p-value", "Masimanimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.O.tr.f)                                        
                                        
write.csv(MasKLPZ.O.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.O.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKLPZ.O.tr.f.status <- merge_samples(MasKLPZ.O.tr.f, MasKLPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.O.tr.f.status <- transform_sample_counts(MasKLPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.O.tr.f.status)                                        
MasKLPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(MasKLPZ.O.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.O.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                 
#KINSHASA vs KONZO HPZ
KinKHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Kinshasa" | KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)                                        
KinKHPZ.O.tr <- transform_sample_counts(KinKHPZ.O, function(x) x / sum(x))                                             
KinKHPZ.O.tr.f <- prune_taxa(f_0.0001, KinKHPZ.O.tr)  

O <- KinKHPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- KinKHPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Kinshasa vs. KHPZ p-value", "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.O.tr.f)                                        
                                        
write.csv(KinKHPZ.O.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.O.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.O.tr.f.status <- merge_samples(KinKHPZ.O.tr.f, KinKHPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.O.tr.f.status <- transform_sample_counts(KinKHPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.O.tr.f.status)                                        
KinKHPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(KinKHPZ.O.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.O.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE) 
                                                 
#MASIMANIMBA AND KONZO HPZ                                        
MasKHPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Masimanimba" | KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)                                        
MasKHPZ.O.tr <- transform_sample_counts(MasKHPZ.O, function(x) x / sum(x))                                             
MasKHPZ.O.tr.f <- prune_taxa(f_0.0001, MasKHPZ.O.tr)  

O <- MasKHPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- MasKHPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "Masimanimba vs. KHPZ p-value", "Masimanimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.O.tr.f)                                        
                                        
write.csv(MasKHPZ.O.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.O.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKHPZ.O.tr.f.status <- merge_samples(MasKHPZ.O.tr.f, MasKHPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.O.tr.f.status <- transform_sample_counts(MasKHPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.O.tr.f.status)                                        
MasKHPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.O.tr.f.status)                                        
                                                                                                
write.csv(t(MasKHPZ.O.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.O.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)                                                                                                  
                                       
#CONTROL (UNAFFECTED)
Control.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.O)                                        
Control.O.tr <- transform_sample_counts(Control.O, function(x) x / sum(x))                                             
Control.O.tr.f <- prune_taxa(f_0.0001, Control.O.tr)  

O <- Control.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- Control.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "ULPZ vs. UHPZ p-value", "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Control_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.O.tr.f.0.05 <- prune_taxa(ls_0.05,Control.O.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.O.tr.f.0.01 <- prune_taxa(ls_0.01,Control.O.tr.f)                                        
                                        
write.csv(Control.O.tr.f.0.05@otu_table, file = "./Control_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.O.tr.f.0.01@otu_table, file = "./Control_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Control.O.tr.f.status <- merge_samples(Control.O.tr.f, Control.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.O.tr.f.status <- transform_sample_counts(Control.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.O.tr.f.status)                                        
Control.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.O.tr.f.status)                                        
                                                                                                
write.csv(t(Control.O.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.O.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                     
#DISEASE (KONZO)
Disease.O <- prune_samples(KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)                                        
Disease.O.tr <- transform_sample_counts(Disease.O, function(x) x / sum(x))                                             
Disease.O.tr.f <- prune_taxa(f_0.0001, Disease.O.tr)  

O <- Disease.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- Disease.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "KLPZ vs. KHPZ p-value", "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Disease_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "Disease_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "Disease_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#Disease.O.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#Disease.O.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.O.tr.f)                                        
                                        
#write.csv(Disease.O.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.O.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Disease.O.tr.f.status <- merge_samples(Disease.O.tr.f, Disease.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Disease.O.tr.f.status <- transform_sample_counts(Disease.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#Disease.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.O.tr.f.status)                                        
#Disease.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.O.tr.f.status)                                        
                                                                                                
#write.csv(t(Disease.O.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.O.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                                     
#LPZ (ULPZ vs. KLPZ)
LPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.O@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.O)                                        
LPZ.O.tr <- transform_sample_counts(LPZ.O, function(x) x / sum(x))                                             
LPZ.O.tr.f <- prune_taxa(f_0.0001, LPZ.O.tr)  

O <- LPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- LPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "ULPZ vs. KLPZ p-value", "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "LPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "LPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "LPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#LPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#LPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.O.tr.f)                                        
                                        
#write.csv(LPZ.O.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.O.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
LPZ.O.tr.f.status <- merge_samples(LPZ.O.tr.f, LPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
LPZ.O.tr.f.status <- transform_sample_counts(LPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.O.tr.f.status)                                        
#LPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.O.tr.f.status)                                        
                                                                                                
#write.csv(t(LPZ.O.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.O.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)
                                             
#HPZ (UHPZ vs. KHPZ)
HPZ.O <- prune_samples(KonzoData.O@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData.O@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.O)                                        
HPZ.O.tr <- transform_sample_counts(HPZ.O, function(x) x / sum(x))                                             
HPZ.O.tr.f <- prune_taxa(f_0.0001, HPZ.O.tr)  

O <- HPZ.O.tr.f
                                               
O.tr_META <- as.data.frame(O@sam_data)
O.tr_OTU <- as.data.frame(t(O@otu_table))
O.tr.DF <- cbind(O.tr_OTU, O.tr_META$Status)

colnames(O.tr.DF)[colnames(O.tr.DF)=="O.tr_META$Status"] <- "Status"
for (i in 1:nrow(O.tr.DF))
  {O.tr.DF[i,]$Status <- HPZ.O.tr.f@sam_data[rownames(O.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(O.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Order", "UHPZ vs. KHPZ p-value", "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(O.tr.DF)-1))
{
  wt <- wilcox.test(O.tr.DF[,i] ~O.tr.DF$Status, data = O.tr.DF)
  WT[i,1] = colnames(O.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "HPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "HPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "HPZ_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#HPZ.O.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.O.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#HPZ.O.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.O.tr.f)                                        
                                        
#write.csv(HPZ.O.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.O.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Order_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
HPZ.O.tr.f.status <- merge_samples(HPZ.O.tr.f, HPZ.O.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
HPZ.O.tr.f.status <- transform_sample_counts(HPZ.O.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.O.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.O.tr.f.status)                                        
#HPZ.O.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.O.tr.f.status)                                        
                                                                                                
#write.csv(t(HPZ.O.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.O.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Order_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                                                                                                                                                                                                                                               
                       
MWW_order <- merge(MWW_order,WT,by="Bacteria Order", sort = FALSE)  
write.csv(MWW_order, file = "Kinshasa_Konzo3_Bacteria_Order_f_0.0001_ByStatus_WilcoxTest_BH.csv")                                              
                                             
#Bacteria Family
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Family")
x <- read.csv("Kinshasa_Konzo3_Family_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                             
#KINSHASA AND MASIMANIMBA
KinMas.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa" | KonzoData.F@sam_data$Status == "Masimanimba", KonzoData.F)                                        
KinMas.F.tr <- transform_sample_counts(KinMas.F, function(x) x / sum(x))                                             
KinMas.F.tr.f <- prune_taxa(f_0.0001, KinMas.F.tr)  

F <- KinMas.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- KinMas.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Kinshasa vs. Masi-manimba p-value", "Kinshasa vs. Masi-manimba p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinMas_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.F.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.F.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.F.tr.f)                                        
                                        
write.csv(KinMas.F.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.F.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.F.tr.f.status <- merge_samples(KinMas.F.tr.f, KinMas.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.F.tr.f.status <- transform_sample_counts(KinMas.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.F.tr.f.status)                                        
KinMas.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.F.tr.f.status)                                        
                                                                                               
write.csv(t(KinMas.F.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.F.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                                                       
MWW_family <- WT
                                                                                               
#KINSHASA AND UNAFFECTED LPZ
KinULPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa" | KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.F)                                        
KinULPZ.F.tr <- transform_sample_counts(KinULPZ.F, function(x) x / sum(x))                                             
KinULPZ.F.tr.f <- prune_taxa(f_0.0001, KinULPZ.F.tr)  

F <- KinULPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- KinULPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Kinshasa vs. ULPZ p-value", "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinULPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinULPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.F.tr.f)                                        
                                        
write.csv(KinULPZ.F.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinULPZ.F.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.F.tr.f.status <- merge_samples(KinULPZ.F.tr.f, KinULPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.F.tr.f.status <- transform_sample_counts(KinULPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinULPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.F.tr.f.status)                                        
#KinULPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(KinULPZ.F.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinULPZ.F.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
                                        
#MASIMANIMBA AND UNAFFECTED LPZ                                        
MasULPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba" | KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.F)                                        
MasULPZ.F.tr <- transform_sample_counts(MasULPZ.F, function(x) x / sum(x))                                             
MasULPZ.F.tr.f <- prune_taxa(f_0.0001, MasULPZ.F.tr)  

F <- MasULPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- MasULPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Masimanimba vs. ULPZ p-value", "Masimanimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.F.tr.f)                                        
                                        
write.csv(MasULPZ.F.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.F.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasULPZ.F.tr.f.status <- merge_samples(MasULPZ.F.tr.f, MasULPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.F.tr.f.status <- transform_sample_counts(MasULPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.F.tr.f.status)                                        
MasULPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(MasULPZ.F.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.F.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)
                                                 
#KINSHASA vs UNAFFECTED HPZ
KinUHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa" | KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.F)                                        
KinUHPZ.F.tr <- transform_sample_counts(KinUHPZ.F, function(x) x / sum(x))                                             
KinUHPZ.F.tr.f <- prune_taxa(f_0.0001, KinUHPZ.F.tr)  

F <- KinUHPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- KinUHPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Kinshasa vs. UHPZ p-value", "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.F.tr.f)                                        
                                        
write.csv(KinUHPZ.F.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.F.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.F.tr.f.status <- merge_samples(KinUHPZ.F.tr.f, KinUHPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.F.tr.f.status <- transform_sample_counts(KinUHPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.F.tr.f.status)                                        
KinUHPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(KinUHPZ.F.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.F.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
                                        
#MASIMANIMBA AND UNAFFECTED HPZ                                        
MasUHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba" | KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.F)                                        
MasUHPZ.F.tr <- transform_sample_counts(MasUHPZ.F, function(x) x / sum(x))                                             
MasUHPZ.F.tr.f <- prune_taxa(f_0.0001, MasUHPZ.F.tr)  

F <- MasUHPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- MasUHPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Masimanimba vs. UHPZ p-value", "Masimanimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "MasUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "MasUHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#MasUHPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#MasUHPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.F.tr.f)                                        
                                        
#write.csv(MasUHPZ.F.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(MasUHPZ.F.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasUHPZ.F.tr.f.status <- merge_samples(MasUHPZ.F.tr.f, MasUHPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.F.tr.f.status <- transform_sample_counts(MasUHPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#MasUHPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.F.tr.f.status)                                        
#MasUHPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.F.tr.f.status)                                        
                                                                                               
#write.csv(t(MasUHPZ.F.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(MasUHPZ.F.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)
                                              
#KINSHASA AND KONZO LPZ
KinKLPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa" | KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.F)                                        
KinKLPZ.F.tr <- transform_sample_counts(KinKLPZ.F, function(x) x / sum(x))                                             
KinKLPZ.F.tr.f <- prune_taxa(f_0.0001, KinKLPZ.F.tr)  

F <- KinKLPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- KinKLPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Kinshasa vs. KLPZ p-value", "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "KinKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKLPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#KinKLPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.F.tr.f)                                        
                                        
write.csv(KinKLPZ.F.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(KinKLPZ.F.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.F.tr.f.status <- merge_samples(KinKLPZ.F.tr.f, KinKLPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.F.tr.f.status <- transform_sample_counts(KinKLPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKLPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.F.tr.f.status)                                        
#KinKLPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(KinKLPZ.F.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(KinKLPZ.F.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
                                        
#MASIMANIMBA AND KONZO LPZ                                        
MasKLPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba" | KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.F)                                        
MasKLPZ.F.tr <- transform_sample_counts(MasKLPZ.F, function(x) x / sum(x))                                             
MasKLPZ.F.tr.f <- prune_taxa(f_0.0001, MasKLPZ.F.tr)  

F <- MasKLPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- MasKLPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Masimanimba vs. KLPZ p-value", "Masimanimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.F.tr.f)                                        
                                        
write.csv(MasKLPZ.F.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.F.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKLPZ.F.tr.f.status <- merge_samples(MasKLPZ.F.tr.f, MasKLPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.F.tr.f.status <- transform_sample_counts(MasKLPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.F.tr.f.status)                                        
MasKLPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(MasKLPZ.F.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.F.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)                                     

#KINSHASA vs KONZO HPZ
KinKHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Kinshasa" | KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)                                        
KinKHPZ.F.tr <- transform_sample_counts(KinKHPZ.F, function(x) x / sum(x))                                             
KinKHPZ.F.tr.f <- prune_taxa(f_0.0001, KinKHPZ.F.tr)  

F <- KinKHPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- KinKHPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Kinshasa vs. KHPZ p-value", "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "KinKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.F.tr.f)                                        
                                        
write.csv(KinKHPZ.F.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.F.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.F.tr.f.status <- merge_samples(KinKHPZ.F.tr.f, KinKHPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.F.tr.f.status <- transform_sample_counts(KinKHPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.F.tr.f.status)                                        
KinKHPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(KinKHPZ.F.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.F.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
                                        
#MASIMANIMBA AND KONZO HPZ                                        
MasKHPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Masimanimba" | KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)                                        
MasKHPZ.F.tr <- transform_sample_counts(MasKHPZ.F, function(x) x / sum(x))                                             
MasKHPZ.F.tr.f <- prune_taxa(f_0.0001, MasKHPZ.F.tr)  

F <- MasKHPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- MasKHPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "Masimanimba vs. KHPZ p-value", "Masimanimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "MasKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.F.tr.f)                                        
                                        
write.csv(MasKHPZ.F.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.F.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKHPZ.F.tr.f.status <- merge_samples(MasKHPZ.F.tr.f, MasKHPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.F.tr.f.status <- transform_sample_counts(MasKHPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.F.tr.f.status)                                        
MasKHPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.F.tr.f.status)                                        
                                                                                               
write.csv(t(MasKHPZ.F.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.F.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)                                                                               
                                       
#CONTROL (UNAFFECTED)
Control.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.F)                                        
Control.F.tr <- transform_sample_counts(Control.F, function(x) x / sum(x))                                             
Control.F.tr.f <- prune_taxa(f_0.0001, Control.F.tr)  

F <- Control.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- Control.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "ULPZ vs. UHPZ p-value", "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Control_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.F.tr.f.0.05 <- prune_taxa(ls_0.05,Control.F.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.F.tr.f.0.01 <- prune_taxa(ls_0.01,Control.F.tr.f)                                        
                                        
write.csv(Control.F.tr.f.0.05@otu_table, file = "./Control_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.F.tr.f.0.01@otu_table, file = "./Control_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Control.F.tr.f.status <- merge_samples(Control.F.tr.f, Control.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.F.tr.f.status <- transform_sample_counts(Control.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.F.tr.f.status)                                        
Control.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.F.tr.f.status)                                        
                                                                                               
write.csv(t(Control.F.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.F.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
                                                 
#DISEASE (KONZO)                                                                     
Disease.F <- prune_samples(KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)                                        
Disease.F.tr <- transform_sample_counts(Disease.F, function(x) x / sum(x))                                             
Disease.F.tr.f <- prune_taxa(f_0.0001, Disease.F.tr)  

F <- Disease.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- Disease.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "KLPZ vs. KHPZ p-value", "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "Disease_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "Disease_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "Disease_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#Disease.F.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#Disease.F.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.F.tr.f)                                        
                                        
#write.csv(Disease.F.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.F.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Disease.F.tr.f.status <- merge_samples(Disease.F.tr.f, Disease.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Disease.F.tr.f.status <- transform_sample_counts(Disease.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#Disease.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.F.tr.f.status)                                        
#Disease.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.F.tr.f.status)                                        
                                                                                               
#write.csv(t(Disease.F.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.F.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)       
                                                 
#LPZ (ULPZ vs. KLPZ)
LPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.F@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.F)                                        
LPZ.F.tr <- transform_sample_counts(LPZ.F, function(x) x / sum(x))                                             
LPZ.F.tr.f <- prune_taxa(f_0.0001, LPZ.F.tr)  

F <- LPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- LPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "ULPZ vs. KLPZ p-value", "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "LPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "LPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "LPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#LPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#LPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.F.tr.f)                                        
                                        
#write.csv(LPZ.F.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.F.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
LPZ.F.tr.f.status <- merge_samples(LPZ.F.tr.f, LPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
LPZ.F.tr.f.status <- transform_sample_counts(LPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.F.tr.f.status)                                        
#LPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.F.tr.f.status)                                        
                                                                                               
#write.csv(t(LPZ.F.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.F.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                 
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE) 
                                                 
#HPZ (UHPZ vs. KHPZ)
HPZ.F <- prune_samples(KonzoData.F@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData.F@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.F)                                        
HPZ.F.tr <- transform_sample_counts(HPZ.F, function(x) x / sum(x))                                             
HPZ.F.tr.f <- prune_taxa(f_0.0001, HPZ.F.tr)  

F <- HPZ.F.tr.f
                                               
F.tr_META <- as.data.frame(F@sam_data)
F.tr_OTU <- as.data.frame(t(F@otu_table))
F.tr.DF <- cbind(F.tr_OTU, F.tr_META$Status)

colnames(F.tr.DF)[colnames(F.tr.DF)=="F.tr_META$Status"] <- "Status"
for (i in 1:nrow(F.tr.DF))
  {F.tr.DF[i,]$Status <- HPZ.F.tr.f@sam_data[rownames(F.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(F.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Family", "UHPZ vs. KHPZ p-value", "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(F.tr.DF)-1))
{
  wt <- wilcox.test(F.tr.DF[,i] ~F.tr.DF$Status, data = F.tr.DF)
  WT[i,1] = colnames(F.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
                                       
WT[,3] <- p.adjust(WT[,2], method = "BH")   
write.csv(WT, file = "HPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                      
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
#write.csv(WT.05, file = "HPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
#write.csv(WT.01, file = "HPZ_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

#ls_0.05 <- WT.05[,1]
#HPZ.F.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.F.tr.f)                                        
#ls_0.01 <- WT.01[,1] 
#HPZ.F.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.F.tr.f)                                        
                                        
#write.csv(HPZ.F.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.F.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Family_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
HPZ.F.tr.f.status <- merge_samples(HPZ.F.tr.f, HPZ.F.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
HPZ.F.tr.f.status <- transform_sample_counts(HPZ.F.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.F.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.F.tr.f.status)                                        
#HPZ.F.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.F.tr.f.status)                                        
                                                                                               
#write.csv(t(HPZ.F.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.F.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Family_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv") 
                                                                                                                                               
MWW_family <- merge(MWW_family,WT,by="Bacteria Family", sort = FALSE)  
write.csv(MWW_family, file = "Kinshasa_Konzo3_Bacteria_Family_f_0.0001_ByStatus_WilcoxTest_BH.csv")                                              

                                                                                             
#Bacteria Genus
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")
x <- read.csv("Kinshasa_Konzo3_Genus_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)

#MWW
                                             
#Kin Mas                                             
KinMas.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Masimanimba"),  KonzoData.G)
KinMas.G.tr <-  transform_sample_counts(KinMas.G, function(x) x / sum(x))
KinMas.G.tr.f <- prune_taxa(f_0.0001, KinMas.G.tr)   
                                        
G <- KinMas.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- KinMas.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Kinshasa vs. Masi-manimba p-value", "Kinshasa vs. Masi-manimba p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinMas_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.G.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.G.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.G.tr.f)                                        
                                        
write.csv(KinMas.G.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.G.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.G.tr.f.status <- merge_samples(KinMas.G.tr.f, KinMas.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.G.tr.f.status <- transform_sample_counts(KinMas.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.G.tr.f.status)                                        
KinMas.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.G.tr.f.status)                                        
                                                
                                                
write.csv(t(KinMas.G.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.G.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

                                                                               
MWW_genus <- WT
                                        
                                        
                                        
#KINSHASA AND UNAFFECTED LPZ
                                             
KinULPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
KinULPZ.G.tr <-  transform_sample_counts(KinULPZ.G, function(x) x / sum(x))
KinULPZ.G.tr.f <- prune_taxa(f_0.0001, KinULPZ.G.tr)  
                                         
#MWW 
                                               
G <- KinULPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- KinULPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Kinshasa vs. ULPZ p-value", "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                           
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinULPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinULPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.G.tr.f)                                        
                                        
write.csv(KinULPZ.G.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinULPZ.G.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.G.tr.f.status <- merge_samples(KinULPZ.G.tr.f, KinULPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.G.tr.f.status <- transform_sample_counts(KinULPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinULPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.G.tr.f.status)                                        
KinULPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(KinULPZ.G.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinULPZ.G.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                         
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)  
                                         
#MASIMANIMBA AND UNAFFECTED LPZ                                        
MasULPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
MasULPZ.G.tr <- transform_sample_counts(MasULPZ.G, function(x) x / sum(x)) 
MasULPZ.G.tr.f <- prune_taxa(f_0.0001, MasULPZ.G.tr)  

#MWW 
                                               
G <- MasULPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- MasULPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Masimanimba vs. ULPZ p-value", "Masimanimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.G.tr.f)                                        
                                        
write.csv(MasULPZ.G.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.G.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasULPZ.G.tr.f.status <- merge_samples(MasULPZ.G.tr.f, MasULPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.G.tr.f.status <- transform_sample_counts(MasULPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.G.tr.f.status)                                        
MasULPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(MasULPZ.G.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.G.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")
                                        
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                        
#KINSHASA AND UNAFFECTED HPZ               
KinUHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
KinUHPZ.G.tr <-  transform_sample_counts(KinUHPZ.G, function(x) x / sum(x))
KinUHPZ.G.tr.f <- prune_taxa(f_0.0001, KinUHPZ.G.tr)  

#MWW 
                                               
G <- KinUHPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- KinUHPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Kinshasa vs. UHPZ p-value", "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.G.tr.f)                                        
                                        
write.csv(KinUHPZ.G.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.G.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.G.tr.f.status <- merge_samples(KinUHPZ.G.tr.f, KinUHPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.G.tr.f.status <- transform_sample_counts(KinUHPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.G.tr.f.status)                                        
KinUHPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(KinUHPZ.G.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.G.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                         
                                         
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                             
                                             
#MASIMANIMBA AND UNAFFECTED HPZ                                        
MasUHPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
MasUHPZ.G.tr <- transform_sample_counts(MasUHPZ.G, function(x) x / sum(x)) 
MasUHPZ.G.tr.f <- prune_taxa(f_0.0001, MasUHPZ.G.tr)  
                                        
#MWW 
                                               
G <- MasUHPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- MasUHPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
                                       
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Masi-manimba vs. UHPZ p-value", "Masi-manimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasUHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasUHPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasUHPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.G.tr.f)                                        
                                        
write.csv(MasUHPZ.G.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasUHPZ.G.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasUHPZ.G.tr.f.status <- merge_samples(MasUHPZ.G.tr.f, MasUHPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.G.tr.f.status <- transform_sample_counts(MasUHPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasUHPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.G.tr.f.status)                                        
MasUHPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(MasUHPZ.G.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasUHPZ.G.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")
                                        
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                             
#KINSHASA and Konzo LPZ
                                             
KinKLPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
KinKLPZ.G.tr <-  transform_sample_counts(KinKLPZ.G, function(x) x / sum(x))
KinKLPZ.G.tr.f <- prune_taxa(f_0.0001, KinKLPZ.G.tr)  
                                         
#MWW 
                                               
G <- KinKLPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- KinKLPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Kinshasa vs. KLPZ p-value", "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKLPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKLPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.G.tr.f)                                        
                                        
write.csv(KinKLPZ.G.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKLPZ.G.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.G.tr.f.status <- merge_samples(KinKLPZ.G.tr.f, KinKLPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.G.tr.f.status <- transform_sample_counts(KinKLPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKLPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.G.tr.f.status)                                        
KinKLPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(KinKLPZ.G.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKLPZ.G.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                                                            
                                         
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)  
                                         
#MASIMANIMBA AND Konzo LPZ                                        
MasKLPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
MasKLPZ.G.tr <- transform_sample_counts(MasKLPZ.G, function(x) x / sum(x)) 
MasKLPZ.G.tr.f <- prune_taxa(f_0.0001, MasKLPZ.G.tr)  

#MWW 
                                               
G <- MasKLPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- MasKLPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Masimanimba vs. KLPZ p-value", "Masimanimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.G.tr.f)                                        
                                        
write.csv(MasKLPZ.G.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.G.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKLPZ.G.tr.f.status <- merge_samples(MasKLPZ.G.tr.f, MasKLPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.G.tr.f.status <- transform_sample_counts(MasKLPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.G.tr.f.status)                                        
MasKLPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(MasKLPZ.G.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.G.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")
                                                                          
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                        
#KINSHASA AND Konzo HPZ               
KinKHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
KinKHPZ.G.tr <-  transform_sample_counts(KinKHPZ.G, function(x) x / sum(x))
KinKHPZ.G.tr.f <- prune_taxa(f_0.0001, KinKHPZ.G.tr)  

#MWW 
                                               
G <- KinKHPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- KinKHPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Kinshasa vs. KHPZ p-value", "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.G.tr.f)                                        
                                        
write.csv(KinKHPZ.G.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.G.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.G.tr.f.status <- merge_samples(KinKHPZ.G.tr.f, KinKHPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.G.tr.f.status <- transform_sample_counts(KinKHPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.G.tr.f.status)                                        
KinKHPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(KinKHPZ.G.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.G.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                         
                                         
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                             
                                             
#MASIMANIMBA AND Konzo HPZ                                        
MasKHPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
MasKHPZ.G.tr <- transform_sample_counts(MasKHPZ.G, function(x) x / sum(x)) 
MasKHPZ.G.tr.f <- prune_taxa(f_0.0001, MasKHPZ.G.tr)  
                                        
#MWW 
                                               
G <- MasKHPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- MasKHPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
                                       
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "Masi-manimba vs. KHPZ p-value", "Masi-manimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.G.tr.f)                                        
                                        
write.csv(MasKHPZ.G.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.G.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKHPZ.G.tr.f.status <- merge_samples(MasKHPZ.G.tr.f, MasKHPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.G.tr.f.status <- transform_sample_counts(MasKHPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.G.tr.f.status)                                        
MasKHPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.G.tr.f.status)                                        
                                                
                                                
write.csv(t(MasKHPZ.G.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.G.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus__WilcoxTest_BH_FDR_0.01.csv")
                                          
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)                                      
                                        
#CONTROL (UNAFFECTED)
Control.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
Control.G.tr <- transform_sample_counts(Control.G, function(x) x / sum(x)) 
Control.G.tr.f <- prune_taxa(f_0.0001, Control.G.tr)  
                                                                                   
G <- Control.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Control.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "ULPZ vs. UHPZ p-value", "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "Control_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.G.tr.f.0.05 <- prune_taxa(ls_0.05,Control.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.G.tr.f.0.01 <- prune_taxa(ls_0.01,Control.G.tr.f)                                        

write.csv(Control.G.tr.f.0.05@otu_table, file = "./Control_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.G.tr.f.0.01@otu_table, file = "./Control_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

Control.G.tr.f.status <- merge_samples(Control.G.tr.f, Control.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.G.tr.f.status <- transform_sample_counts(Control.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.G.tr.f.status)                                        
Control.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.G.tr.f.status)                                        

write.csv(t(Control.G.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.G.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
                                              
   
#DISEASE (KONZO)
Disease.G <- prune_samples((KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
Disease.G.tr <- transform_sample_counts(Disease.G, function(x) x / sum(x)) 
Disease.G.tr.f <- prune_taxa(f_0.0001, Disease.G.tr)  

# Konzo LPZ vs. HPZ
G <- Disease.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Disease.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
                                            
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "KLPZ vs. KHPZ p-value", "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "Disease_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Disease_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Disease_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#Disease.G.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#Disease.G.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.G.tr.f)                                        

#write.csv(Disease.G.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.G.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

#Disease.G.tr.f.status <- merge_samples(Disease.G.tr.f, Disease.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
#Disease.G.tr.f.status <- transform_sample_counts(Disease.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#Disease.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.G.tr.f.status)                                        
#Disease.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.G.tr.f.status)                                        

#write.csv(t(Disease.G.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.G.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)

#Low Prevalence Zone (LPZ)

LPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
LPZ.G.tr <- transform_sample_counts(LPZ.G, function(x) x / sum(x)) 
LPZ.G.tr.f <- prune_taxa(f_0.0001, LPZ.G.tr)  

                                                
G <- LPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- LPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "ULPZ vs. KLPZ p-value", "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "LPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "LPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "LPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#LPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#LPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.G.tr.f)                                        

#write.csv(LPZ.G.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.G.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

#LPZ.G.tr.f.status <- merge_samples(LPZ.G.tr.f, LPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
#LPZ.G.tr.f.status <- transform_sample_counts(LPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.G.tr.f.status)                                        
#LPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.G.tr.f.status)                                        

#write.csv(t(LPZ.G.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.G.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                    
                                    
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)

                                                
#High Prevalence Zone (HPZ)
HPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
HPZ.G.tr <- transform_sample_counts(HPZ.G, function(x) x / sum(x)) 
HPZ.G.tr.f <- prune_taxa(f_0.0001, HPZ.G.tr)  

                                               
G <- HPZ.G.tr.f
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- HPZ.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(G.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Genus", "UHPZ vs. KHPZ p-value", "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(G.tr.DF)-1))
{
  wt <- wilcox.test(G.tr.DF[,i] ~G.tr.DF$Status, data = G.tr.DF)
  WT[i,1] = colnames(G.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "HPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                    
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "HPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "HPZ_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#HPZ.G.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.G.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#HPZ.G.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.G.tr.f)                                        

#write.csv(HPZ.G.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.G.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Genus_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        

#HPZ.G.tr.f.status <- merge_samples(HPZ.G.tr.f, HPZ.G.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
#HPZ.G.tr.f.status <- transform_sample_counts(HPZ.G.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.G.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.G.tr.f.status)                                        
#HPZ.G.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.G.tr.f.status)                                        

#write.csv(t(HPZ.G.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.G.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Genus_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                    
                                    
MWW_genus <- merge(MWW_genus,WT,by="Bacteria Genus", sort = FALSE)
write.csv(MWW_genus, file = "Kinshasa_Konzo3_Bacteria_Genus_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                                       
#Bacteria Species
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Species")
x <- read.csv("Kinshasa_Konzo3_Species_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                             
#KINSHASA AND MASIMANIMBA
KinMas.S <-  prune_samples(KonzoData.S@sam_data$Status == "Kinshasa" | KonzoData.S@sam_data$Status == "Masimanimba", KonzoData.S)
KinMas.S.tr <-  transform_sample_counts(KinMas.S, function(x) x / sum(x))
KinMas.S.tr.f <- prune_taxa(f_0.0001, KinMas.S.tr)  
                                        
#MWW 
#Make a data frame with the rel abund values for the samples being compared, and add Status (Group) as a column for testing                                               
S <- KinMas.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- KinMas.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Kinshasa vs. Masi-manimba p-value", "Kinshasa vs. Masi-manimba p-value adjusted")
for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinMas_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinMas_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinMas_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinMas.S.tr.f.0.05 <- prune_taxa(ls_0.05,KinMas.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinMas.S.tr.f.0.01 <- prune_taxa(ls_0.01,KinMas.S.tr.f)                                        
                                        
write.csv(KinMas.S.tr.f.0.05@otu_table, file = "./KinMas_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinMas.S.tr.f.0.01@otu_table, file = "./KinMas_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinMas.S.tr.f.status <- merge_samples(KinMas.S.tr.f, KinMas.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinMas.S.tr.f.status <- transform_sample_counts(KinMas.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinMas.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinMas.S.tr.f.status)                                        
KinMas.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinMas.S.tr.f.status)                                        
                                                                                              
write.csv(t(KinMas.S.tr.f.status.0.05@otu_table), file = "./KinMas_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinMas.S.tr.f.status.0.01@otu_table), file = "./KinMas_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                      
                                        
MWW_species <- WT                                              

#KINSHASA AND UNAFFECTED LPZ
                                             
KinULPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa" | KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.S)
KinULPZ.S.tr <-  transform_sample_counts(KinULPZ.S, function(x) x / sum(x))
KinULPZ.S.tr.f <- prune_taxa(f_0.0001, KinULPZ.S.tr)  

                                               
S <- KinULPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- KinULPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Kinshasa vs. ULPZ p-value", "Kinshasa vs. ULPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinULPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,KinULPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinULPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,KinULPZ.S.tr.f)                                        
                                        
write.csv(KinULPZ.S.tr.f.0.05@otu_table, file = "./KinULPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinULPZ.S.tr.f.0.01@otu_table, file = "./KinULPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinULPZ.S.tr.f.status <- merge_samples(KinULPZ.S.tr.f, KinULPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinULPZ.S.tr.f.status <- transform_sample_counts(KinULPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinULPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinULPZ.S.tr.f.status)                                        
KinULPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinULPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(KinULPZ.S.tr.f.status.0.05@otu_table), file = "./KinULPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinULPZ.S.tr.f.status.0.01@otu_table), file = "./KinULPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                            
                                         
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)
                                              

#MASIMANIMBA AND UNAFFECTED LPZ                                        
MasULPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba" | KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData.S)
MasULPZ.S.tr <- transform_sample_counts(MasULPZ.S, function(x) x / sum(x)) 
MasULPZ.S.tr.f <- prune_taxa(f_0.0001, MasULPZ.S.tr)  

#MWW 
                                               
S <- MasULPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- MasULPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Masi-manimba vs. ULPZ p-value",  "Masi-manimba vs. ULPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasULPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasULPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,MasULPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasULPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,MasULPZ.S.tr.f)                                        
                                        
write.csv(MasULPZ.S.tr.f.0.05@otu_table, file = "./MasULPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasULPZ.S.tr.f.0.01@otu_table, file = "./MasULPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasULPZ.S.tr.f.status <- merge_samples(MasULPZ.S.tr.f, MasULPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasULPZ.S.tr.f.status <- transform_sample_counts(MasULPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasULPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasULPZ.S.tr.f.status)                                        
MasULPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasULPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(MasULPZ.S.tr.f.status.0.05@otu_table), file = "./MasULPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasULPZ.S.tr.f.status.0.01@otu_table), file = "./MasULPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                 
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

#KINSHASA vs UHPZ (Kin vs. CI)
KinUHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa" | KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.S)
KinUHPZ.S.tr <-  transform_sample_counts(KinUHPZ.S, function(x) x / sum(x))
KinUHPZ.S.tr.f <- prune_taxa(f_0.0001, KinUHPZ.S.tr)  


#MWW                                       
S <- KinUHPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- KinUHPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Kinshasa vs. UHPZ p-value", "Kinshasa vs. UHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinUHPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,KinUHPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinUHPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,KinUHPZ.S.tr.f)                                        
                                        
write.csv(KinUHPZ.S.tr.f.0.05@otu_table, file = "./KinUHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinUHPZ.S.tr.f.0.01@otu_table, file = "./KinUHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinUHPZ.S.tr.f.status <- merge_samples(KinUHPZ.S.tr.f, KinUHPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinUHPZ.S.tr.f.status <- transform_sample_counts(KinUHPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinUHPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinUHPZ.S.tr.f.status)                                        
KinUHPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinUHPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(KinUHPZ.S.tr.f.status.0.05@otu_table), file = "./KinUHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinUHPZ.S.tr.f.status.0.01@otu_table), file = "./KinUHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")  
                                                 
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

                                         
#MASIMANIMBA vs. UHPZ
MasUHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba" | KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.S)
MasUHPZ.S.tr <- transform_sample_counts(MasUHPZ.S, function(x) x / sum(x)) 
MasUHPZ.S.tr.f <- prune_taxa(f_0.0001, MasUHPZ.S.tr)  
                                       
#MWW 
                                               
S <- MasUHPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- MasUHPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    

WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Masi-manimba vs. UHPZ p-value",  "Masi-manimba vs. UHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasUHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasUHPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,MasUHPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasUHPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,MasUHPZ.S.tr.f)                                        
                                        
write.csv(MasUHPZ.S.tr.f.0.05@otu_table, file = "./MasUHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasUHPZ.S.tr.f.0.01@otu_table, file = "./MasUHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasUHPZ.S.tr.f.status <- merge_samples(MasUHPZ.S.tr.f, MasUHPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasUHPZ.S.tr.f.status <- transform_sample_counts(MasUHPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasUHPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasUHPZ.S.tr.f.status)                                        
MasUHPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasUHPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(MasUHPZ.S.tr.f.status.0.05@otu_table), file = "./MasUHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasUHPZ.S.tr.f.status.0.01@otu_table), file = "./MasUHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

#KINSHASA AND KONZO LPZ
                                             
KinKLPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa" | KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.S)
KinKLPZ.S.tr <-  transform_sample_counts(KinKLPZ.S, function(x) x / sum(x))
KinKLPZ.S.tr.f <- prune_taxa(f_0.0001, KinKLPZ.S.tr)  

                                               
S <- KinKLPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- KinKLPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Kinshasa vs. KLPZ p-value", "Kinshasa vs. KLPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKLPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,KinKLPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKLPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,KinKLPZ.S.tr.f)                                        
                                        
write.csv(KinKLPZ.S.tr.f.0.05@otu_table, file = "./KinKLPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKLPZ.S.tr.f.0.01@otu_table, file = "./KinKLPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKLPZ.S.tr.f.status <- merge_samples(KinKLPZ.S.tr.f, KinKLPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKLPZ.S.tr.f.status <- transform_sample_counts(KinKLPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKLPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKLPZ.S.tr.f.status)                                        
KinKLPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKLPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(KinKLPZ.S.tr.f.status.0.05@otu_table), file = "./KinKLPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKLPZ.S.tr.f.status.0.01@otu_table), file = "./KinKLPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                            
                                                                                 
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)
                                              

#MASIMANIMBA AND KONZO LPZ                                        
MasKLPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba" | KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.S)
MasKLPZ.S.tr <- transform_sample_counts(MasKLPZ.S, function(x) x / sum(x)) 
MasKLPZ.S.tr.f <- prune_taxa(f_0.0001, MasKLPZ.S.tr)  

#MWW 
                                               
S <- MasKLPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- MasKLPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Masi-manimba vs. KLPZ p-value",  "Masi-manimba vs. KLPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKLPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKLPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,MasKLPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKLPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,MasKLPZ.S.tr.f)                                        
                                        
write.csv(MasKLPZ.S.tr.f.0.05@otu_table, file = "./MasKLPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKLPZ.S.tr.f.0.01@otu_table, file = "./MasKLPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKLPZ.S.tr.f.status <- merge_samples(MasKLPZ.S.tr.f, MasKLPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKLPZ.S.tr.f.status <- transform_sample_counts(MasKLPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKLPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKLPZ.S.tr.f.status)                                        
MasKLPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKLPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(MasKLPZ.S.tr.f.status.0.05@otu_table), file = "./MasKLPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKLPZ.S.tr.f.status.0.01@otu_table), file = "./MasKLPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                                                    
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

#KINSHASA vs KHPZ
KinKHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Kinshasa" | KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)
KinKHPZ.S.tr <-  transform_sample_counts(KinKHPZ.S, function(x) x / sum(x))
KinKHPZ.S.tr.f <- prune_taxa(f_0.0001, KinKHPZ.S.tr)  


#MWW                                       
S <- KinKHPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- KinKHPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Kinshasa vs. KHPZ p-value", "Kinshasa vs. KHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "KinKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                         
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "KinKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "KinKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
KinKHPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,KinKHPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
KinKHPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,KinKHPZ.S.tr.f)                                        
                                        
write.csv(KinKHPZ.S.tr.f.0.05@otu_table, file = "./KinKHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(KinKHPZ.S.tr.f.0.01@otu_table, file = "./KinKHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
KinKHPZ.S.tr.f.status <- merge_samples(KinKHPZ.S.tr.f, KinKHPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
KinKHPZ.S.tr.f.status <- transform_sample_counts(KinKHPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

KinKHPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,KinKHPZ.S.tr.f.status)                                        
KinKHPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,KinKHPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(KinKHPZ.S.tr.f.status.0.05@otu_table), file = "./KinKHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(KinKHPZ.S.tr.f.status.0.01@otu_table), file = "./KinKHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                            
                                         
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

                                         
#MASIMANIMBA vs. KHPZ
MasKHPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Masimanimba" | KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)
MasKHPZ.S.tr <- transform_sample_counts(MasKHPZ.S, function(x) x / sum(x)) 
MasKHPZ.S.tr.f <- prune_taxa(f_0.0001, MasKHPZ.S.tr)  
                                       
#MWW 
                                               
S <- MasKHPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- MasKHPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    

WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "Masi-manimba vs. KHPZ p-value",  "Masi-manimba vs. KHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "MasKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "MasKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "MasKHPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
MasKHPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,MasKHPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
MasKHPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,MasKHPZ.S.tr.f)                                        
                                        
write.csv(MasKHPZ.S.tr.f.0.05@otu_table, file = "./MasKHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(MasKHPZ.S.tr.f.0.01@otu_table, file = "./MasKHPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
MasKHPZ.S.tr.f.status <- merge_samples(MasKHPZ.S.tr.f, MasKHPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
MasKHPZ.S.tr.f.status <- transform_sample_counts(MasKHPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

MasKHPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,MasKHPZ.S.tr.f.status)                                        
MasKHPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,MasKHPZ.S.tr.f.status)                                        
                                                                                              
write.csv(t(MasKHPZ.S.tr.f.status.0.05@otu_table), file = "./MasKHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(MasKHPZ.S.tr.f.status.0.01@otu_table), file = "./MasKHPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")
                                                                                                                     
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)
                                    
#CONTROL (UNAFFECTED)
Control.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData.S)
Control.S.tr <- transform_sample_counts(Control.S, function(x) x / sum(x)) 
Control.S.tr.f <- prune_taxa(f_0.0001, Control.S.tr)  
                                        
                                               
S <- Control.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- Control.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "ULPZ vs. UHPZ p-value", "ULPZ vs. UHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "Control_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")

WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Control_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Control_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
Control.S.tr.f.0.05 <- prune_taxa(ls_0.05,Control.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
Control.S.tr.f.0.01 <- prune_taxa(ls_0.01,Control.S.tr.f)                                        
                                        
write.csv(Control.S.tr.f.0.05@otu_table, file = "./Control_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
write.csv(Control.S.tr.f.0.01@otu_table, file = "./Control_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Control.S.tr.f.status <- merge_samples(Control.S.tr.f, Control.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Control.S.tr.f.status <- transform_sample_counts(Control.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

Control.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,Control.S.tr.f.status)                                        
Control.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,Control.S.tr.f.status)                                        
                                                                                              
write.csv(t(Control.S.tr.f.status.0.05@otu_table), file = "./Control_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
write.csv(t(Control.S.tr.f.status.0.01@otu_table), file = "./Control_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                             
                                        
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)

                                      
#DISEASE (KONZO)
Disease.S <- prune_samples(KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)                                        
Disease.S.tr <- transform_sample_counts(Disease.S, function(x) x / sum(x))                                             
Disease.S.tr.f <- prune_taxa(f_0.0001, Disease.S.tr)  

S <- Disease.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- Disease.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "KLPZ vs. KHPZ p-value", "KLPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "Disease_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                        
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "Disease_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "Disease_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#Disease.S.tr.f.0.05 <- prune_taxa(ls_0.05,Disease.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#Disease.S.tr.f.0.01 <- prune_taxa(ls_0.01,Disease.S.tr.f)                                        
                                        
#write.csv(Disease.S.tr.f.0.05@otu_table, file = "./Disease_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(Disease.S.tr.f.0.01@otu_table, file = "./Disease_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
Disease.S.tr.f.status <- merge_samples(Disease.S.tr.f, Disease.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
Disease.S.tr.f.status <- transform_sample_counts(Disease.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#Disease.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,Disease.S.tr.f.status)                                        
#Disease.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,Disease.S.tr.f.status)                                        
                                                                                              
#write.csv(t(Disease.S.tr.f.status.0.05@otu_table), file = "./Disease_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(Disease.S.tr.f.status.0.01@otu_table), file = "./Disease_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                             
                                        
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)


#Low Prevalence Zone (LPZ)
LPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData.S@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData.S)                                        
LPZ.S.tr <- transform_sample_counts(LPZ.S, function(x) x / sum(x))                                             
LPZ.S.tr.f <- prune_taxa(f_0.0001, LPZ.S.tr)  
                                             
S <- LPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- LPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
    
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "ULPZ vs. KLPZ p-value", "ULPZ vs. KLPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "LPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                    
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "LPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "LPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#LPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,LPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#LPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,LPZ.S.tr.f)                                        
                                        
#write.csv(LPZ.S.tr.f.0.05@otu_table, file = "./LPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(LPZ.S.tr.f.0.01@otu_table, file = "./LPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
#LPZ.S.tr.f.status <- merge_samples(LPZ.S.tr.f, LPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
#LPZ.S.tr.f.status <- transform_sample_counts(LPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#LPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,LPZ.S.tr.f.status)                                        
#LPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,LPZ.S.tr.f.status)                                        
                                                                                              
#write.csv(t(LPZ.S.tr.f.status.0.05@otu_table), file = "./LPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(LPZ.S.tr.f.status.0.01@otu_table), file = "./LPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                             
                                    
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)
                                                
#High Prevalence Zone (HPZ)
HPZ.S <- prune_samples(KonzoData.S@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData.S@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData.S)                                        
HPZ.S.tr <- transform_sample_counts(HPZ.S, function(x) x / sum(x))                                             
HPZ.S.tr.f <- prune_taxa(f_0.0001, HPZ.S.tr)  

S <- HPZ.S.tr.f
                                               
S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
for (i in 1:nrow(S.tr.DF))
  {S.tr.DF[i,]$Status <- HPZ.S.tr.f@sam_data[rownames(S.tr.DF[i,]),]$Status
  }
                                                
WT <- matrix(nrow = ncol(S.tr_OTU), ncol = 3)
colnames(WT) <- c("Bacteria Species", "UHPZ vs. KHPZ p-value",  "UHPZ vs. KHPZ p-value adjusted")

for (i in 1:(ncol(S.tr.DF)-1))
{
  wt <- wilcox.test(S.tr.DF[,i] ~S.tr.DF$Status, data = S.tr.DF)
  WT[i,1] = colnames(S.tr.DF[i])
  WT[i,2] = as.numeric(wt$p.value)
}
WT[,3] <- p.adjust(WT[,2], method = "BH")                                          
write.csv(WT, file = "HPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                    
WT.05 <- subset(WT, as.numeric(WT[,3]) <= 0.05)
write.csv(WT.05, file = "HPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
WT.01 <- subset(WT, as.numeric(WT[,3]) <= 0.01)
write.csv(WT.01, file = "HPZ_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH_FDR_0.01.csv")

ls_0.05 <- WT.05[,1]
#HPZ.S.tr.f.0.05 <- prune_taxa(ls_0.05,HPZ.S.tr.f)                                        
ls_0.01 <- WT.01[,1] 
#HPZ.S.tr.f.0.01 <- prune_taxa(ls_0.01,HPZ.S.tr.f)                                        
                                        
#write.csv(HPZ.S.tr.f.0.05@otu_table, file = "./HPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")
#write.csv(HPZ.S.tr.f.0.01@otu_table, file = "./HPZ_Bacteria_Species_f_0.0001_RelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                        
                                        
#HPZ.S.tr.f.status <- merge_samples(HPZ.S.tr.f, HPZ.S.tr.f@sam_data$Status) #merge_smaples by default sums the values for otu
#HPZ.S.tr.f.status <- transform_sample_counts(HPZ.S.tr.f.status, function(x) x / 30) #average the sum of relabund in each group

#HPZ.S.tr.f.status.0.05 <- prune_taxa(ls_0.05,HPZ.S.tr.f.status)                                        
#HPZ.S.tr.f.status.0.01 <- prune_taxa(ls_0.01,HPZ.S.tr.f.status)                                        
                                                                                              
#write.csv(t(HPZ.S.tr.f.status.0.05@otu_table), file = "./HPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.05.csv")                                                
#write.csv(t(HPZ.S.tr.f.status.0.01@otu_table), file = "./HPZ_Bacteria_Species_f_0.0001_AvgRelAbund_ByStatus_WilcoxTest_BH_FDR_0.01.csv")                                                                                                              
                                    
MWW_species <- merge(MWW_species,WT,by="Bacteria Species", sort = FALSE)
write.csv(MWW_species, file = "Kinshasa_Konzo3_Bacteria_Species_f_0.0001_ByStatus_WilcoxTest_BH.csv")
                                                                                                            
### Beta Diversity using Bray-Curtis for Bacteria Genus ----------------- 
                                             
#Using phyloseq distance function and bray method, and the sample wise distances are generated. Using the adonis (PERMANOVA) function, stats are performed with 99999 permutations
#the p-values generated here are reported in PCoA Figures for the relavant comparisons  

setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")
x <- read.csv("Kinshasa_Konzo3_Genus_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)  
                         

KonzoData.G.tr.f <- prune_taxa(f_0.0001, KonzoData.G.tr)

#Initially perform a multivariate analysis to determine which variables influence the gut microbiome of the various cohorts                                    
otuD.G <- as.data.frame(t(otu_table(KonzoData.G.f)))
diversity.G <- estimate_richness(KonzoData.G.f)
diversity.G <- cbind(sample_data(KonzoData.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
brayd <- phyloseq::distance(KonzoData.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Geography * diversity.G$Region * diversity.G$Disease * diversity.G$Age * diversity.G$Sex, perm=99999); bdiv_bray

#Specific Comparisons  

Geography.G <- prune_samples((KonzoData.G@sam_data$Status != "Konzo_Low_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_High_Prevalence_Zone"), KonzoData.G)                                              
Geography.G.tr <-  transform_sample_counts(Geography.G, function(x) x / sum(x))
Geography.G.tr.log10 <-  transform_sample_counts(Geography.G.tr, function(x) log10(x))    
Geography.G.f <- prune_taxa(f_0.0001, Geography.G)                                                 
Geography.G.tr.f <- prune_taxa(f_0.0001, Geography.G.tr)


otuD.G <- as.data.frame(t(otu_table(Geography.G.f)))
diversity.G <- estimate_richness(Geography.G.f)
diversity.G <- cbind(sample_data(Geography.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(Geography.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_Geography.txt")
                                                 

KinMas.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Masimanimba"),  KonzoData.G)
KinMas.G.tr <-  transform_sample_counts(KinMas.G, function(x) x / sum(x))
KinMas.G.f <- prune_taxa(f_0.0001, KinMas.G) 
KinMas.G.tr.f <- prune_taxa(f_0.0001, KinMas.G.tr)   


otuD.G <- as.data.frame(t(otu_table(KinMas.G.f)))
diversity.G <- estimate_richness(KinMas.G.f)
diversity.G <- cbind(sample_data(KinMas.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Kinshasa", "Masimanimba"))

brayd <- phyloseq::distance(KinMas.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_KinMas.txt")


KinULPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
KinULPZ.G.tr <-  transform_sample_counts(KinULPZ.G, function(x) x / sum(x))
KinULPZ.G.f <- prune_taxa(f_0.0001, KinULPZ.G)
KinULPZ.G.tr.f <- prune_taxa(f_0.0001, KinULPZ.G.tr)  


otuD.G <- as.data.frame(t(otu_table(KinULPZ.G.f)))
diversity.G <- estimate_richness(KinULPZ.G.f)
diversity.G <- cbind(sample_data(KinULPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Kinshasa", "Unaffected_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(KinULPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_KinULPZ.txt")


MasULPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone"),  KonzoData.G)
MasULPZ.G.tr <- transform_sample_counts(MasULPZ.G, function(x) x / sum(x)) 
MasULPZ.G.f <- prune_taxa(f_0.0001, MasULPZ.G)
MasULPZ.G.tr.f <- prune_taxa(f_0.0001, MasULPZ.G.tr)  

otuD.G <- as.data.frame(t(otu_table(MasULPZ.G.f)))
diversity.G <- estimate_richness(MasULPZ.G.f)
diversity.G <- cbind(sample_data(MasULPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Masimanimba", "Unaffected_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(MasULPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_MasULPZ.txt")

KinUHPZ.G <-  prune_samples((KonzoData.G@sam_data$Status == "Kinshasa") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
KinUHPZ.G.tr <-  transform_sample_counts(KinUHPZ.G, function(x) x / sum(x))
KinUHPZ.G.f <- prune_taxa(f_0.0001, KinUHPZ.G)
KinUHPZ.G.tr.f <- prune_taxa(f_0.0001, KinUHPZ.G.tr)  

otuD.G <- as.data.frame(t(otu_table(KinUHPZ.G.f)))
diversity.G <- estimate_richness(KinUHPZ.G.f)
diversity.G <- cbind(sample_data(KinUHPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Kinshasa", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(KinUHPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_KinUHPZ.txt")



MasUHPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Masimanimba") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
MasUHPZ.G.tr <- transform_sample_counts(MasUHPZ.G, function(x) x / sum(x)) 
MasUHPZ.G.f <- prune_taxa(f_0.0001, MasUHPZ.G)
MasUHPZ.G.tr.f <- prune_taxa(f_0.0001, MasUHPZ.G.tr)  

otuD.G <- as.data.frame(t(otu_table(MasUHPZ.G.f)))
diversity.G <- estimate_richness(MasUHPZ.G.f)
diversity.G <- cbind(sample_data(MasUHPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Masimanimba", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(MasUHPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_MasUHPZ.txt")
  
Control.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone"),  KonzoData.G)
Control.G.tr <- transform_sample_counts(Control.G, function(x) x / sum(x)) 
Control.G.f <- prune_taxa(f_0.0001, Control.G)                                         
Control.G.tr.f <- prune_taxa(f_0.0001, Control.G.tr)                                          
                                        
                                        
otuD.G <- as.data.frame(t(otu_table(Control.G.f)))
diversity.G <- estimate_richness(Control.G.f)
diversity.G <- cbind(sample_data(Control.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(Control.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_Control.txt")

                                        
Disease.G <- prune_samples((KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
Disease.G.tr <- transform_sample_counts(Disease.G, function(x) x / sum(x)) 
Disease.G.f <- prune_taxa(f_0.0001, Disease.G)  
Disease.G.tr.f <- prune_taxa(f_0.0001, Disease.G.tr)                                                                                 
                                        
otuD.G <- as.data.frame(t(otu_table(Disease.G.f)))
diversity.G <- estimate_richness(Disease.G.f)
diversity.G <- cbind(sample_data(Disease.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(Disease.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_Disease.txt")
                                        
LPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_Low_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_Low_Prevalence_Zone"),  KonzoData.G)
LPZ.G.tr <- transform_sample_counts(LPZ.G, function(x) x / sum(x)) 
LPZ.G.f <- prune_taxa(f_0.0001, LPZ.G)  
LPZ.G.tr.f <- prune_taxa(f_0.0001, LPZ.G.tr)                                                                                                                       

otuD.G <- as.data.frame(t(otu_table(LPZ.G.f)))
diversity.G <- estimate_richness(LPZ.G.f)
diversity.G <- cbind(sample_data(LPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(LPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_LPZ.txt")
                                    

HPZ.G <- prune_samples((KonzoData.G@sam_data$Status == "Unaffected_High_Prevalence_Zone") | (KonzoData.G@sam_data$Status == "Konzo_High_Prevalence_Zone"),  KonzoData.G)
HPZ.G.tr <- transform_sample_counts(HPZ.G, function(x) x / sum(x)) 
HPZ.G.f <- prune_taxa(f_0.0001, HPZ.G)                                      
HPZ.G.tr.f <- prune_taxa(f_0.0001, HPZ.G.tr)                                     
                                    
otuD.G <- as.data.frame(t(otu_table(HPZ.G.f)))
diversity.G <- estimate_richness(HPZ.G.f)
diversity.G <- cbind(sample_data(HPZ.G.f),diversity.G) #Might change since cbind can be tricky and not reliable, so always confirm if correctly done
diversity.G$Status <- as.factor(diversity.G$Status)
diversity.G$Status <- factor(diversity.G$Status, levels = c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(HPZ.G.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ diversity.G$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_filtered_bdiv_genus_adonis_HPZ.txt")


#Figure 3: Geography Genus Filtered PCoA

#Genus and Corr PCoA Figure
#Extract eigen values (values of the variance reported) after running ordinate function with bray distance (used by phyloseq in generating PCoA plots)
b <- ordinate(Geography.G.tr.f, method="PCoA", dist="bray")
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

G <- Geography.G.tr.f
G.tr.DF <- as.data.frame(t(G@otu_table))

n = ncol(G.tr.DF)

G.tr.DF["Status"] <- NA
G.tr.DF["Axis.1"] <- NA
G.tr.DF["Axis.2"] <- NA

G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Geography.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
  }

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Axis.1 <- b.DF[rownames(G.tr.DF[i,]),1]
  }

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Axis.2 <- b.DF[rownames(G.tr.DF[i,]),2]
  }

Cor1 <- matrix(nrow = n,  ncol = 3)

colnames(Cor1) <- c("Genus vs. Axis.1", "spearman cor", "p-value")
for (i in 1:n)
{
  cor <- cor.test(G.tr.DF[,i], G.tr.DF$Axis.1, method=c("spearman"))
  Cor1[i,1] = colnames(G.tr.DF[i])
  Cor1[i,2] = as.numeric(cor$estimate)
  Cor1[i,3] = as.numeric(cor$p.value)
}
write.csv(Cor1, file = "Geography_Genus_f_0.0001_Axis.1_Correlation.csv")
Cor2 <- matrix(nrow = n,  ncol = 3)

colnames(Cor2) <- c("Genus vs. Axis.2", "spearman cor", "p-value")
for (i in 1:n)
{
  cor <- cor.test(G.tr.DF[,i], G.tr.DF$Axis.2, method=c("spearman"))
  Cor2[i,1] = colnames(G.tr.DF[i])
  Cor2[i,2] = as.numeric(cor$estimate)
  Cor2[i,3] = as.numeric(cor$p.value)
}
write.csv(Cor2, file = "Geography_Genus_f_0.0001_Axis.2_Correlation.csv")                                        
                                        
#Plot most correlated with Axis 1 and Axis 2                                             
#Correlation Plot
#Axis 1 Prevotella
#Axis 2 Faecalibacterium

a1 <- ggplot(G.tr.DF, aes(x = Axis.1, y = Prevotella)) +
    geom_point(aes(color = factor(Status)), size = 1, stroke = 0, shape = 16) + theme_classic() + ylab("Prev.") + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 7))
a1 <- a1 + scale_color_manual(labels = SL, values = geography_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.y = element_text(size = 7, face = "italic"), axis.text.y = element_text(size = 7))  
a1 <- a1 + geom_smooth(method=lm, color = "black", size = 0.5) + theme(plot.margin=unit(c(0.15,0.15,0.25,0.15), "lines")) + scale_y_continuous(breaks = seq(-0.1, 0.4, by = 0.2))
#a1 <- a1 + stat_cor(method = "spearman", size = 5) 

a2 <- ggplot(G.tr.DF, aes(x = Faecalibacterium, y = Axis.2)) +
    geom_point(aes(color = factor(Status)), size = 1, stroke = 0, shape = 16) + theme_classic() + xlab("Faec.") + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 7))
a2 <- a2 + scale_color_manual(labels = SL, values = geography_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.x = element_text(size = 7, face = "italic"), axis.text.x = element_text(size = 7))  
a2 <- a2 + scale_y_continuous(position = "right") + scale_x_continuous(position = "top", breaks = seq(0.1, 0.3, by = 0.2))
a2 <- a2 + geom_smooth(method=lm, color = "black", size = 0.5) + theme(plot.margin=unit(c(0.15,0.125,0.15,0.15), "lines"))

p1 = plot_ordination(Geography.G.tr.f, ordinate(Geography.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
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
PGBt <- PGBt + annotate("text", x = -0.46, y = -0.42, label = expression(paste("p = 1x",10^-5)), size = 2.5) #1e-05
PGBt <- ggarrange(PGBt,labels = c("A"),font.label = list(size = 7))

PGB <- PGB + scale_x_continuous(position = "top") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines")) + theme(legend.position="none")


G <- arrangeGrob(PGBt, a1,                               # bar plot spaning two columns
             a2, l,                               # box plot and scatter plot
             ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1,1,1,3), c(1,1,1,3), c(1,1,1,3), c(2, 2, 2, 4)))

tiff(filename = "Overall_Geography_Genus_Filtered_PCoA_Corr.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
ggarrange(as_ggplot(G))
dev.off()

#KinMas
p1 = plot_ordination(KinMas.G.tr.f, ordinate(KinMas.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PKMB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinmas_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

PKMBt <- PKMB + stat_ellipse(type = "t") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
PKMBt <- PKMBt + theme(legend.position="none")
PKMBt <- PKMBt + annotate("text", x = 0.4, y = -0.43, label = expression(paste("p = 2x",10^-5)), size = 2) #2e-05
PKMBt <- ggarrange(PKMBt,labels = c("B"),font.label = list(size = 7))
                                    
#KinULPZ
p1 = plot_ordination(KinULPZ.G.tr.f, ordinate(KinULPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PKUB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinulpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

PKUBt <- PKUB + stat_ellipse(type = "t") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
PKUBt <- PKUBt + theme(legend.position="none")
PKUBt <- PKUBt + annotate("text", x = 0.37, y = -0.43, label = expression(paste("p = 0.00139")), size = 2)#0.00139
PKUBt <- ggarrange(PKUBt,labels = c("C"),font.label = list(size = 7))
                                    
#MasULPZ
p1 = plot_ordination(MasULPZ.G.tr.f, ordinate(MasULPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PMUB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = masulpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

PMUBt <- PMUB + stat_ellipse(type = "t") + scale_x_continuous(position = "top") + scale_y_continuous(position = "right") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
PMUBt <- PMUBt + theme(legend.position="none")
PMUBt <- PMUBt + annotate("text", x = 0.38, y = -0.27, label = expression(paste("p = 3x",10^-5)), size = 2) #3e-05
PMUBt <- ggarrange(PMUBt,labels = c("E"),font.label = list(size = 7))

#KinUHPZ
p1 = plot_ordination(KinUHPZ.G.tr.f, ordinate(KinUHPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PKUHB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinuhpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

PKUHBt <- PKUHB + stat_ellipse(type = "t") +  scale_y_continuous(position = "right")+ theme(plot.margin=unit(c(0.15,0.15,0.15,0.6), "lines"))
PKUHBt <- PKUHBt + theme(legend.position="none")
PKUHBt <- PKUHBt + annotate("text", x = 0.38, y = -0.43, label = expression(paste("p = 1x",10^-5)), size = 2) #1e-05
PKUHBt <- ggarrange(PKUHBt,labels = c("D"),font.label = list(size = 7))

#MasUHPZ
p1 = plot_ordination(MasUHPZ.G.tr.f, ordinate(MasUHPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
PMUHB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = masuhpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

PMUHBt <- PMUHB + stat_ellipse(type = "t") + scale_y_continuous(position = "right") + theme(plot.margin=unit(c(0.6,0.15,0.15,0.15), "lines"))
PMUHBt <- PMUHBt + theme(legend.position="none")
PMUHBt <- PMUHBt + annotate("text", x = 0.28, y = -0.33, label = expression(paste("p = 0.00321")), size = 2) #0.00321
PMUHBt <- ggarrange(PMUHBt,labels = c("F"),font.label = list(size = 7))

#EF <- ggarrange(PMUBt, PMUHBt, ncol = 1, nrow = 2, labels = c("E", "F"))
#BCD <- ggarrange(PKMBt, PKUBt, PKUHBt, ncol = 3, nrow = 1, align = "h")
#Geo_Test <- ggarrange(G, EF, BCD, ncol = 2, nrow = 2, widths = c(3,1), heights = c(3,1))
#V1
Geo <- arrangeGrob(PGBt, a1,                              
             a2, l, PMUBt, PMUHBt, PKMBt, PKUBt, PKUHBt,                             
             ncol = 6, nrow = 6,
             layout_matrix = rbind(c(1,1,1,3,5,5), c(1,1,1,3,5,5), c(1,1,1,3,6,6), c(2,2,2,4,6,6), c(7,7,8,8,9,9), c(7,7,8,8,9,9)))

tiff(filename = "Geography_Genus_Filtered_PCoA_Corr.tiff", width = 5.5, height = 5.5, units = "in", res = 600)
ggarrange(as_ggplot(Geo))
dev.off()
                                        

                                             
###Figure 5: Kahemba Prevalence Zones Genus PCoA------------------------

b <- ordinate(Control.G.tr.f, method="PCoA", dist="bray")
b.DF <- as.data.frame(b$vectors)

G <- Control.G.tr.f
G.tr.DF <- as.data.frame(t(G@otu_table))

n = ncol(G.tr.DF)

G.tr.DF["Status"] <- NA
G.tr.DF["Axis.1"] <- NA
G.tr.DF["Axis.2"] <- NA

G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Control.G.tr@sam_data[rownames(G.tr.DF[i,]),]$Status
  }

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Axis.1 <- b.DF[rownames(G.tr.DF[i,]),1]
  }

for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Axis.2 <- b.DF[rownames(G.tr.DF[i,]),2]
  }

Cor1 <- matrix(nrow = n,  ncol = 3)

colnames(Cor1) <- c("Genus vs. Axis.1", "spearman cor", "p-value")
for (i in 1:n)
{
  cor <- cor.test(G.tr.DF[,i], G.tr.DF$Axis.1, method=c("spearman"))
  Cor1[i,1] = colnames(G.tr.DF[i])
  Cor1[i,2] = as.numeric(cor$estimate)
  Cor1[i,3] = as.numeric(cor$p.value)
}
write.csv(Cor1, file = "Control_Genus_f_0.0001_Axis.1_Correlation.csv")
Cor2 <- matrix(nrow = n,  ncol = 3)

colnames(Cor2) <- c("Genus vs. Axis.2", "spearman cor", "p-value")
for (i in 1:n)
{
  cor <- cor.test(G.tr.DF[,i], G.tr.DF$Axis.2, method=c("spearman"))
  Cor2[i,1] = colnames(G.tr.DF[i])
  Cor2[i,2] = as.numeric(cor$estimate)
  Cor2[i,3] = as.numeric(cor$p.value)
}
write.csv(Cor2, file = "Control_Genus_f_0.0001_Axis.2_Correlation.csv")
                                    
#Correlation Plots
#Faecalibacterium Axis 1
#Prevotella Axis 2

a1 <- ggplot(G.tr.DF, aes(x = Axis.1, y = Faecalibacterium)) +
    geom_point(aes(color = factor(Status)), size = 0.7, stroke = 0, shape = 16) + theme_classic() + ylab("Faec.") + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 5))
a1 <- a1 + scale_color_manual(labels = SSSL, values = control_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.y = element_text(size = 6, face = "italic"), axis.text.y = element_text(size = 5))  
a1 <- a1 + geom_smooth(method=lm, color = "black", size = 0.35) + theme(plot.margin=unit(c(0.05,0.05,0.05,0.05), "lines")) + scale_y_continuous(breaks = seq(-0.1, 0.3, by = 0.2))
#a1 <- a1 + stat_cor(method = "spearman", size = 5) 

a2 <- ggplot(G.tr.DF, aes(x = Prevotella, y = Axis.2)) +
    geom_point(aes(color = factor(Status)), size = 0.7, stroke = 0, shape = 16) + theme_classic() + xlab("Prev.") + theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 5))
a2 <- a2 + scale_color_manual(labels = SSSL, values = control_color) + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme (axis.title.x = element_text(size = 6, face = "italic"), axis.text.x = element_text(size = 5))  
a2 <- a2 + scale_y_continuous(position = "right") + scale_x_continuous(position = "top", breaks = seq(0.1, 0.4, by = 0.2))
a2 <- a2 + geom_smooth(method=lm, color = "black", size = 0.35) + theme(plot.margin=unit(c(0.05,0.05,0.05,0.05), "lines"))

                                    
p1 = plot_ordination(Control.G.tr.f, ordinate(Control.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 0.85, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]

PCB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = control_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.2, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 6), axis.title.x = element_text(size = 6), axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5))

P <- PCB + guides(colour = guide_legend(override.aes = list(size=1)))

l <- get_legend(P)
l <- as_ggplot(l)
l <- l + theme(plot.margin=unit(c(-1,0,0,-1), "lines"))


PCBt <- PCB + stat_ellipse(type = "t") + scale_x_continuous(position = "top", breaks = seq(-0.1, 0.5, by = 0.5) ) + theme(plot.margin=unit(c(0.05,0.05,0.05,0.05), "lines"))


PCBt <- PCBt + theme(legend.position="none")
PCBt <- PCBt + annotate("text", x = 0.39, y = -0.24, label = expression(paste("p = 0.00057")), size = 2) #0.00057

C <- arrangeGrob(PCBt, a1,                               # bar plot spaning two columns
             a2, l,                               # box plot and scatter plot
             ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1,1,1,3), c(1,1,1,3), c(1,1,1,3), c(2, 2, 2, 4)))
as_ggplot(C)

#Prevotella and Faecalibacterium for ULPZ vs. UHPZ #relative abundance

G <- Control.G.tr.f

G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)
G.tr.DF <- cbind(G.tr.DF, G.tr_META$Geography)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Geography"] <- "Geography"

for (i in nrow(G.tr.DF))
{G.tr.DF[i,]$Status <- Control.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Status
}
G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Unaffected_Low_Prevalence_Zone" ,"Unaffected_High_Prevalence_Zone"))

for (i in nrow(G.tr.DF))
{G.tr.DF[i,]$Geography <- Control.G.tr.f@sam_data[rownames(G.tr.DF[i,]),]$Geography
}
G.tr.DF$Geography <- factor(G.tr.DF$Geography, levels = c("Low_Prevalence_Zone", "High_Prevalence_Zone"))

#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

my_comparisons <- list(c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

#remove outlier.shape = NA and add outlier.size if you don't want jitter
p <- ggplot(G.tr.DF,aes(x = Status,y = Prevotella)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab(expression(paste("rel. abund. of ", italic("Prevotella")))) + stat_boxplot(geom ='errorbar')
p <- p + geom_jitter(position=position_jitter(0.2), size = 0.35)
p <- p + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = control_color) + theme(plot.title = element_blank(), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 6), axis.title.x = element_blank())
p <- p + stat_summary(fun.y=mean, geom="point", shape=23, size=1, color="black", fill="white")
#remove outlier.shape = NA and add outlier.size if you don't want jitter
f <- ggplot(G.tr.DF,aes(x = Status,y = Faecalibacterium)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab(expression(paste("rel. abund. of ", italic("Faecalibacterium")))) + stat_boxplot(geom ='errorbar')
f <- f + geom_jitter(position=position_jitter(0.2), size = 0.35)
f <- f + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = control_color) + theme(plot.title = element_blank(), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 6), axis.title.x = element_blank())
f <- f + stat_summary(fun.y=mean, geom="point", shape=23, size=1, color="black", fill="white")
                                    
#Prevotella and Faecalibacterium for ULPZ vs. UHPZ #CLR
p2 <- ggplot(Control.G.CLR.DF,aes(x = Status,y = Prevotella)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab(expression(paste("median CLR value of ", italic("Prevotella")))) + stat_boxplot(geom ='errorbar')
p2 <- p2 + geom_jitter(position=position_jitter(0.2), size = 0.35)
p2 <- p2 + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = control_color) + theme(plot.title = element_blank(), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 6), axis.title.x = element_blank())
p2 <- p2 + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
p2 <- p2 + stat_summary(fun.y=mean, geom="point", shape=23, size=1, color="black", fill="white")
                                    
f2 <- ggplot(Control.G.CLR.DF,aes(x = Status,y = Faecalibacterium)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab(expression(paste("median CLR value of ", italic("Faecalibacterium")))) + stat_boxplot(geom ='errorbar')
f2 <- f2 + geom_jitter(position=position_jitter(0.2), size = 0.35)
f2 <- f2 + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = control_color) + theme(plot.title = element_blank(), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 4), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 6), axis.title.x = element_blank())
f2 <- f2 + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
f2 <- f2 + stat_summary(fun.y=mean, geom="point", shape=23, size=1, color="black", fill="white")
                                    
                                    
pf <- ggarrange(p,f, p2, f2, ncol = 4, nrow = 1, align = "hv")                                  
fp <- ggarrange (f, p, f2, p2, ncol = 4, nrow = 1, align = "hv")
                                   
part1 <- ggarrange(as_ggplot(C), pf, labels = c("A","B"), font.label = list(size = 7), ncol = 2, nrow = 1, widths = c(2.5, 3.5)) 
                                    
part1_v2 <- ggarrange(as_ggplot(C),fp, labels = c("A","B"), font.label = list(size = 7), ncol = 2, nrow = 1, widths = c(2.5, 3.5))
                                    
#KLPZ vs. KHPZ

p1 = plot_ordination(Disease.G.tr.f, ordinate(Disease.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
PKB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = disease_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom",legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

PKB <- PKB + guides(colour = guide_legend(override.aes = list(size=1)))

PKBt <- PKB + stat_ellipse(type = "t", show.legend = FALSE) 

PKBt <- PKBt + annotate("text", x = 0.4, y = -0.34, label = "p = 0.01744", size = 2) #0.01744

PKBt                                    
                                    
#ULPZ vs. KLPZ
p1 = plot_ordination(LPZ.G.tr.f, ordinate(LPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
PNIB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = lpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

PNIBt <- PNIB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
PNIBt <- PNIBt + annotate("text", x = 0.42, y = -0.56, label = expression(paste("p = 0.9105")), size = 2) #0.9105
                                    
#UHPZ vs. KHPZ
p1 = plot_ordination(HPZ.G.tr.f, ordinate(HPZ.G.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
PIB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = hpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))


PIBt <- PIB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
PIBt <- PIBt + annotate("text", x = 0.5, y = -0.4, label = expression(paste("p = 0.569")), size = 2) #0.569

                                    
part2 <- ggarrange(PKBt, PNIBt, PIBt, labels = c("C","D", "E"), ncol = 3, nrow = 1, font.label = list(size = 7), widths = c(2, 2, 2) )                                  

tiff(filename = "Kahemba_Genus_Filtered_PCoA_CombFig.tiff", width = 7, height = 4.5, units = "in", res = 600)
ggarrange(part1, part2, ncol = 1, nrow = 2, heights = c(2.5,2))
dev.off()  

tiff(filename = "Kahemba_Genus_Filtered_PCoA_CombFig_V2.tiff", width = 7, height = 4.5, units = "in", res = 600)
ggarrange(part1_v2, part2, ncol = 1, nrow = 2, heights = c(2.5,2))
dev.off()  

                                   
                                    
##### Supplemental Figures                                   
                                    
#Supplementary Figure 1
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken")                                             
diet <- read.csv("./Kinshasa_Konzo3_Diet_3.csv")

diet$Status <- factor(diet$Status, levels = c("Kinshasa", "Masimanimba", "Low_Prevalence_Zone", "High_Prevalence_Zone"))
diet$Food <- factor(diet$Food, levels = c("Maize/Rice", "Other Cereals", "Bread/Wheat", "Tubers and Other Roots", "Pulses", "Meat", "Oil Crops", "Milk/Dairy Products", "Vegetables", "Other Greens", "Fruits", "Sweets", "Red Palm Oil", "Vegetable Oil"))

empty_bar <- 2
to_add <- data.frame( matrix(NA, empty_bar*nlevels(diet$Status), ncol(diet)) )
colnames(to_add) <- colnames(diet)
to_add$Status <- rep(levels(diet$Status), each=empty_bar)
diet <- rbind(diet, to_add)
diet <- diet %>% arrange(Status)
diet$id <- seq(1, nrow(diet))

label_diet <- diet
number_of_bar <- nrow(label_diet)
angle <- 90 - 360 * (label_diet$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_diet$hjust <- ifelse( angle < -90, 1, 0)
label_diet$angle <- ifelse(angle < -90, angle+180, angle)

base_diet <- diet %>% 
  dplyr::group_by(Status) %>% 
  dplyr::summarize(start=min(id), end=max(id) - empty_bar) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_diet <- base_diet
grid_diet$end <- grid_diet$end[ c(nrow(grid_diet), 1:nrow(grid_diet)-1)] + 1
grid_diet$start <- grid_diet$start - 1
grid_diet <- grid_diet[-1,]


p <- ggplot(diet, aes(x=as.factor(id), y=Frequency, fill=Food)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  geom_bar(aes(x=as.factor(id), y=Frequency, fill=Food), colour = "black", stat="identity", width = 1) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_diet, aes(x = end, y = 7, xend = start, yend = 7), colour = "black", alpha=1, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_diet, aes(x = end, y = 5, xend = start, yend = 5), colour = "black", alpha=1, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_diet, aes(x = end, y = 3, xend = start, yend = 3), colour = "black", alpha=1, size=0.2 , inherit.aes = FALSE ) +
  geom_segment(data=grid_diet, aes(x = end, y = 1, xend = start, yend = 1), colour = "black", alpha=1, size=0.2 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(diet$id),4), y = c(1, 3, 5, 7), label = c("1d", "3d", "5d", "7d") , color="black", size=3 , angle=0, hjust=1) +
  
  ylim(-10,NA) +
  theme_minimal() +
  theme(
    legend.position = c(0.5,0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm"), 
    legend.text = element_text(size = 5),
    legend.key.size = unit(.25, "cm"),
    legend.title = element_blank()
  ) + 
  guides(fill=guide_legend(ncol=2,byrow=TRUE)) +
  coord_polar() +
  geom_segment(data=base_diet, aes(x = start, y = -1, xend = end, yend = -1), colour = "black", size=0.3 , inherit.aes = FALSE ) + geom_text(data=base_diet, aes(x = title, y = -2.5, label=c("Kinshasa", "Masi-Manimba","LPZ", "HPZ")), hjust=c(0.5,0.5, 0.5,0.5), angle=c(-45, 45, -45, 45), colour = "black", size=2.5, inherit.aes = FALSE)

p

tiff(filename = "KinshasaKonzo3_Diet.tiff", width = 5.5, height = 5.5, units = "in", res = 600)
p
dev.off()                                     
                                             
## Supplementary Figure 2:

#Phylum
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Phylum")
                                     
top_P <- read.csv("Kinshasa_Konzo3_Phylum_Top4.csv", row.names = 1, colClasses = "character")
top_P <- unlist(top_P)

KonzoData.P.tr.status.top = prune_taxa(top_P, KonzoData.P.tr.status)

p <- plot_bar(KonzoData.P.tr.status.top, "Sample", "Abundance", fill = 'phylum')
p$data$Sample <- factor(p$data$Sample, levels = c("Konzo_High_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_Low_Prevalence_Zone", "Masimanimba", "Kinshasa"))

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
                                    
#Class
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Class")

top_C <- read.csv("Kinshasa_Konzo3_Class_Top5.csv", row.names = 1, colClasses = "character")
top_C <- unlist(top_C)

KonzoData.C.tr.status.top = prune_taxa(top_C, KonzoData.C.tr.status)

p <- plot_bar(KonzoData.C.tr.status.top, "Sample", "Abundance", fill = 'class')
p$data$Sample <- factor(p$data$Sample, levels = c("Konzo_High_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_Low_Prevalence_Zone", "Masimanimba", "Kinshasa"))
#p <- p + geom_bar(stat = "identity") + coord_flip()
p <- p + labs(x = element_blank(), y = "Relative Abundance") +  scale_fill_discrete(name = "Class")
top_class_plot <- p + theme(legend.position="bottom") + theme(legend.key.size = unit(.4, "cm"))


top_class_plot <- top_class_plot + 
  theme(legend.position="right", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-8,0,-8,-8)) + 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.9)) +
  scale_x_discrete(labels= SSSL)+
  theme(plot.title = element_blank(), legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 7), legend.title = element_text(size = 7)) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme(axis.text.y = element_text(angle = 0, size = 7), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 0,vjust=1, hjust=0.5, size = 7), axis.title.x = element_text(size = 7), legend.title = element_text(size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_text(size=7))
top_class_plot <- top_class_plot + geom_bar(stat = "identity") + coord_flip() + theme(axis.title.y = element_blank())
top_class_plot

#Order
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Order")

top_O <- read.csv("Kinshasa_Konzo3_Order_Top5.csv", row.names = 1, colClasses = "character")
top_O <- unlist(top_O)

KonzoData.O.tr.status.top = prune_taxa(top_O, KonzoData.O.tr.status)

p <- plot_bar(KonzoData.O.tr.status.top, "Sample", "Abundance", fill = 'order')
p$data$Sample <- factor(p$data$Sample, levels = c("Konzo_High_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_Low_Prevalence_Zone", "Masimanimba", "Kinshasa"))

p <- p + labs(x = element_blank(), y = "Relative Abundance") +  scale_fill_discrete(name = "Order")
top_order_plot <- p + theme(legend.position="bottom") + theme(legend.key.size = unit(.4, "cm"))

top_order_plot <- top_order_plot + 
  theme(legend.position="right", legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-8,0,-8,-8))+ 
  scale_y_continuous(expand = c(0,0), limits = c(0,0.9)) +
  scale_x_discrete(labels= SSSL)+
  theme(plot.title = element_blank(), legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 7), legend.title = element_text(size = 7)) + 
  guides(fill=guide_legend(ncol=1,byrow=TRUE)) + 
  theme(axis.text.y = element_text(angle = 0, size = 7), axis.ticks.y = element_blank(), axis.text.x = element_text(angle = 0,vjust=1, hjust=0.5, size = 7), axis.title.x = element_text(size = 7), legend.title = element_text(size = 7)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), strip.text.x = element_text(size=7))
top_order_plot <- top_order_plot + geom_bar(stat = "identity") + coord_flip() + theme(axis.title.y = element_blank())
top_order_plot
                                    

#Family
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Family")
                                     
top_F <- read.csv("Kinshasa_Konzo3_Family_Top5.csv", row.names = 1, colClasses = "character")
top_F <- unlist(top_F)

KonzoData.F.tr.status.top = prune_taxa(top_F, KonzoData.F.tr.status)

p <- plot_bar(KonzoData.F.tr.status.top, "Sample", "Abundance", fill = 'family')
p$data$Sample <- factor(p$data$Sample, levels = c("Konzo_High_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_Low_Prevalence_Zone", "Masimanimba", "Kinshasa"))

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

                                    
#Box Plot for Ratio of Prevotella / Bacteroides (Genus)
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")
                                    
#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone")) 


temp <- as.data.frame(t(otu_table(KonzoData.G.tr)))
which( colnames(temp)=="Prevotella" )
BP <- temp[,400:401, drop = FALSE]

BP$Status <- "Mis"
BP$Ratio <- 0
BP$Status <- factor(BP$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone","Konzo_Low_Prevalence_Zone","Unaffected_High_Prevalence_Zone","Konzo_High_Prevalence_Zone"))

for (i in 1:nrow(BP))
{BP[i,]$Status <- KonzoData.G.tr@sam_data[rownames(BP[i,]),]$Status
}

for (i in 1:nrow(BP))
{BP[i,]$Ratio <- BP[i,1]/BP[i,2]
}

bp1 <- ggboxplot(BP, x = "Status", y = "Ratio" , fill = "Status", xlab = "Samples", ylab = "rel. abund. of Prevotella/Bacteroides", title = "", outlier.size = 1)
bp2 <- bp1 + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "t.test", size = 2)
bp3 <- bp2 + theme(legend.position="none") + theme(axis.title.y = element_text(size = 7), axis.title.x = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) + stat_boxplot(geom ="errorbar") + scale_fill_manual(values = konzo_color) + scale_x_discrete(labels = SSSL)
bp3
#tiff(filename = "KinshasaKonzo3_Genus_RelAbund_Prev_over_Bact.tiff", width = 2.5, height = 3, units = "in", res = 600)
#bp3
#dev.off()

#Average Bray-Curtis for Genus with Relative Abundance (NOT READ COUNTS) for Each Group (Intra)

brayd <- phyloseq::distance(KonzoData.G.tr.f, method="bray")
brayd.DF <- as.data.frame(as.matrix(brayd))

KIN <- matrix(nrow = 30, ncol = 30)
colnames(KIN) <- colnames(brayd.DF[1:30])
rownames(KIN) <- rownames(brayd.DF[1:30,])
for (i in 1:30) {for (j in 1:30) {KIN[i,j] <- brayd.DF[i,j]}}
write.csv(KIN, file = "Kinshasa_RelAbundBray_Genus_filtered.csv")

KI <- matrix(nrow = 30, ncol = 30)
colnames(KI) <- colnames(brayd.DF[31:60])
rownames(KI) <- rownames(brayd.DF[31:60,])
for (i in 1:30) {for (j in 1:30) {KI[i,j] <- brayd.DF[i+30,j+30]}}
write.csv(KI, file = "KHPZ_RelAbundBray_Genus_filtered.csv")

CI <- matrix(nrow = 30, ncol = 30)
colnames(CI) <- colnames(brayd.DF[61:90])
rownames(CI) <- rownames(brayd.DF[61:90,])
for (i in 1:30) {for (j in 1:30) {CI[i,j] <- brayd.DF[i+60,j+60]}}
write.csv(CI, file = "UHPZ_RelAbundBray_Genus_filtered.csv")

KNI <- matrix(nrow = 30, ncol = 30)
colnames(KNI) <- colnames(brayd.DF[91:120])
rownames(KNI) <- rownames(brayd.DF[91:120,])
for (i in 1:30) {for (j in 1:30) {KNI[i,j] <- brayd.DF[i+90,j+90]}}
write.csv(KNI, file = "KLPZ_RelAbundBray_Genus_filtered.csv")

CNI <- matrix(nrow = 30, ncol = 30)
colnames(CNI) <- colnames(brayd.DF[121:150])
rownames(CNI) <- rownames(brayd.DF[121:150,])
for (i in 1:30) {for (j in 1:30) {CNI[i,j] <- brayd.DF[i+120,j+120]}}
write.csv(CNI, file = "ULPZ_RelAbundBray_Genus_filtered.csv")

MAS <- matrix(nrow = 30, ncol = 30)
colnames(MAS) <- colnames(brayd.DF[151:180])
rownames(MAS) <- rownames(brayd.DF[151:180,])
for (i in 1:30) {for (j in 1:30) {MAS[i,j] <- brayd.DF[i+150,j+150]}}
write.csv(MAS, file = "Masimanimba_RelAbundBray_Genus_filtered.csv")

bray_avg <- matrix(nrow = 180,  ncol = 3)

colnames(bray_avg) <- c("Samples", "BrayAvg", "Status")

for (i in 1:30) 
{
  bray_avg[i,1] <- rownames(KIN)[i]
  bray_avg[i,2] <- mean(KIN[i,])
  bray_avg[i,3] <- "Kinshasa"
}

for (i in 1:30) 
{
  bray_avg[i+30,1] <- rownames(MAS)[i]
  bray_avg[i+30,2] <- mean(MAS[i,])
  bray_avg[i+30,3] <- "Masimanimba"
}

for (i in 1:30) 
{
  bray_avg[i+60,1] <- rownames(CNI)[i]
  bray_avg[i+60,2] <- mean(CNI[i,])
  bray_avg[i+60,3] <- "Unaffected_Low_Prevalence_Zone"
}

for (i in 1:30) 
{
  bray_avg[i+90,1] <- rownames(KNI)[i]
  bray_avg[i+90,2] <- mean(KNI[i,])
  bray_avg[i+90,3] <- "Konzo_Low_Prevalence_Zone"
}

for (i in 1:30) 
{
  bray_avg[i+120,1] <- rownames(CI)[i]
  bray_avg[i+120,2] <- mean(CI[i,])
  bray_avg[i+120,3] <- "Unaffected_High_Prevalence_Zone"
}

for (i in 1:30) 
{
  bray_avg[i+150,1] <- rownames(KI)[i]
  bray_avg[i+150,2] <- mean(KI[i,])
  bray_avg[i+150,3] <- "Konzo_High_Prevalence_Zone"
}

write.csv(bray_avg, file = "KinshasaControl_Konzo3_RelAbundBray_Averages_Genus_filtered_PerSamples.csv")

bray_avg <- read.csv("./KinshasaControl_Konzo3_RelAbundBray_Averages_Genus_filtered_PerSamples.csv")

bray_avg$Status <- factor(bray_avg$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))


b <- ggplot(bray_avg, aes(factor(Status), BrayAvg)) + geom_boxplot(aes(fill = factor(Status)), outlier.size = 1) + labs(x = element_blank(), y = "Average Bray-Curtis") + theme(axis.text.x = element_blank()) + theme_classic()
b2 <- b + theme(legend.position="none")
b3 <- b2 + theme(axis.title.y = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7)) + stat_boxplot(geom ="errorbar")
b4 <- b3 + scale_fill_manual(values = konzo_color) + scale_x_discrete(labels = SSSL)
#tiff("KinshasaControl_Konzo3_Bacteria_Genus_AvgRelAbundBrayDistancePerSamples.tiff", width = 2, height = 2.5, units = "in", res = 600)
#b4
#dev.off()

#Supplemental Fig 2: Class, Order, Prev/Bact and Intra Bray
tiff(filename = "Konzo1Konzo3_P_C_O_F_PrevOverBact_IntraBrayFiltered.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(top_phylum_plot, top_class_plot, top_order_plot, top_family_plot, bp3, b4, labels = c("A","B","C", "D", "E", "F"), ncol = 2, nrow = 3, font.label = list(size = 7))
dev.off()
                               
                                                                        
##Supplementary Fig 4: Kin vs. Mas vs. ULPZ
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken/Bacteria/Bacteria_Genus")

#Random Forest Box Plots
s <- plot_spacer() + theme_minimal()
#spacer plot is added so the percent breakdown for top ten taxa can be added using gimp

##----IMPORTANT CHANGE----------- Geography (Kin vs. Mas vs. ULPZ only)
                                     
my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_Low_Prevalence_Zone")) 
g_color <- c("royalblue1",   "springgreen3", "turquoise3")

#Specifically for this Figure
                                     
Geography.G <- prune_samples((KonzoData.G@sam_data$Status != "Unaffected_High_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_Low_Prevalence_Zone") & (KonzoData.G@sam_data$Status != "Konzo_High_Prevalence_Zone"), KonzoData.G)                                              
Geography.G.tr <-  transform_sample_counts(Geography.G, function(x) x / sum(x))
Geography.G.tr.log10 <-  transform_sample_counts(Geography.G.tr, function(x) log10(x))

G <- Geography.G.tr.log10
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Geography.G.tr.log10@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone"))

#Kinshasa vs. ULPZ vs. Mas
                                     
#Actinomyces                                     
#Clostridioides
#Leuconostoc                                     
#Megasphaera                                     
#Mageeibacillus 
#Selenomonas                                                  
#Kordia                                                                       
#Arachidicoccus
#Cohnella
#Aneurinibacillus
                                                                                                                                                                  
                                    
g1 <- ggboxplot(G.tr.DF, x = "Status", y = "Actinomyces", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Actinomyces"))) + theme(axis.text.x = element_blank()) + theme_classic()
g1 <- g1 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g1 <- g1 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g1 <- g1 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g2 <- ggboxplot(G.tr.DF, x = "Status", y = "Clostridioides", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Clostridioides"))) + theme(axis.text.x = element_blank()) + theme_classic()
g2 <- g2 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g2 <- g2 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g2 <- g2 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g3 <- ggboxplot(G.tr.DF, x = "Status", y = "Leuconostoc", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Leuconostoc"))) + theme(axis.text.x = element_blank()) + theme_classic()
g3 <- g3 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g3 <- g3 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g3 <- g3 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g4 <- ggboxplot(G.tr.DF, x = "Status", y = "Megasphaera", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Megasphaera"))) + theme(axis.text.x = element_blank()) + theme_classic()
g4 <- g4 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g4 <- g4 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g4 <- g4 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g5 <- ggboxplot(G.tr.DF, x = "Status", y = "Mageeibacillus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Mageeibacillus"))) + theme(axis.text.x = element_blank()) + theme_classic()
g5 <- g5 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g5 <- g5 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g5 <- g5 + guides(fill=guide_legend(ncol=3,byrow=TRUE))
                                                 
g6 <- ggboxplot(G.tr.DF, x = "Status", y = "Selenomonas", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Selenomonas"))) + theme(axis.text.x = element_blank()) + theme_classic()
g6 <- g6 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g6 <- g6 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g6 <- g6 + guides(fill=guide_legend(ncol=3,byrow=TRUE))                                                 

g7 <- ggboxplot(G.tr.DF, x = "Status", y = "Kordia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Kordia"))) + theme(axis.text.x = element_blank()) + theme_classic()
g7 <- g7 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g7 <- g7 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g7 <- g7 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g8 <- ggboxplot(G.tr.DF, x = "Status", y = "Arachidicoccus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Arachidicoccus"))) + theme(axis.text.x = element_blank()) + theme_classic()
g8 <- g8 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g8 <- g8 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g8 <- g8 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g9 <- ggboxplot(G.tr.DF, x = "Status", y = "Cohnella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Cohnella"))) + theme(axis.text.x = element_blank()) + theme_classic()
g9 <- g9 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g9 <- g9 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g9 <- g9 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g10 <- ggboxplot(G.tr.DF, x = "Status", y = "Aneurinibacillus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Aneurinibacillus"))) + theme(axis.text.x = element_blank()) + theme_classic()
g10 <- g10 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g10 <- g10 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g10 <- g10 + guides(fill=guide_legend(ncol=3,byrow=TRUE))                                                                                                 

tiff(filename = "KinshasaKonzo3_Genus_RF_Kin_Geo_Boxplot.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(s, s, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, ncol = 4, nrow = 3, labels = c("", "", "A","B","C","D","E","F","G","H","I","J"), font.label = list(size = 6), common.legend = TRUE, legend = "bottom")
dev.off()    
                
## Supp 5 Mas vs. (Kin and ULPZ)
                                                 
                                                 
#Masimanimba vs. All
                                     
#Phoenicibacter                                     
#Tolumonas                                                                         
#Rothia 
#Faecalibacterium
#Collinsella                                     
#Salmonella
#Megasphaera  
#Faecalitalea                                                  
#Actinomyces                                                                        
#Gordonibacter                                                 
                                                 
g11 <- ggboxplot(G.tr.DF, x = "Status", y = "Phoenicibacter", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Phoenicibacter"))) + theme(axis.text.x = element_blank()) + theme_classic()
g11 <- g11 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g11 <- g11 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g11 <- g11 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g12 <- ggboxplot(G.tr.DF, x = "Status", y = "Tolumonas", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Tolumonas"))) + theme(axis.text.x = element_blank()) + theme_classic()
g12 <- g12 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g12 <- g12 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g12 <- g12 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g13 <- ggboxplot(G.tr.DF, x = "Status", y = "Rothia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Rothia"))) + theme(axis.text.x = element_blank()) + theme_classic()
g13 <- g13 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g13 <- g13 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g13 <- g13 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g14 <- ggboxplot(G.tr.DF, x = "Status", y = "Faecalibacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Faecalibacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
g14 <- g14 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g14 <- g14 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g14 <- g14 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g15 <- ggboxplot(G.tr.DF, x = "Status", y = "Collinsella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Collinsella"))) + theme(axis.text.x = element_blank()) + theme_classic()
g15 <- g15 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g15 <- g15 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g15 <- g15 + guides(fill=guide_legend(ncol=3,byrow=TRUE))                                                 
                                                 
g16 <- ggboxplot(G.tr.DF, x = "Status", y = "Salmonella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Salmonella"))) + theme(axis.text.x = element_blank()) + theme_classic()
g16 <- g16 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g16 <- g16 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g16 <- g16 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g17 <- ggboxplot(G.tr.DF, x = "Status", y = "Megasphaera", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Megasphaera"))) + theme(axis.text.x = element_blank()) + theme_classic()
g17 <- g17 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g17 <- g17 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g17 <- g17 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g18 <- ggboxplot(G.tr.DF, x = "Status", y = "Faecalitalea", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Faecalitalea"))) + theme(axis.text.x = element_blank()) + theme_classic()
g18 <- g18 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g18 <- g18 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g18 <- g18 + guides(fill=guide_legend(ncol=3,byrow=TRUE))
                                                 
g19 <- ggboxplot(G.tr.DF, x = "Status", y = "Actinomyces", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Actinomyces"))) + theme(axis.text.x = element_blank()) + theme_classic()
g19 <- g19 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g19 <- g19 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g19 <- g19 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g20 <- ggboxplot(G.tr.DF, x = "Status", y = "Gordonibacter", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Gordonibacter"))) + theme(axis.text.x = element_blank()) + theme_classic()
g20 <- g20 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g20 <- g20 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g20 <- g20 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

tiff(filename = "KinshasaKonzo3_Genus_RF_Mas_Geo_Boxplot.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(s,s,g11, g12, g13, g14, g15, g16, g17, g18, g19, g20, ncol = 4, nrow = 3, labels = c("", "", "A","B","C","D","E","F","G","H","I","J"), font.label = list(size = 6), common.legend = TRUE, legend = "bottom")
dev.off()                                                                                                 
                                                               
## Supp 6: ULPZ vs. (Kin and Mas)
                                                 
#ULPZ vs. All   
                                     
#Denitrobacterium    
#Gemmantimonas                                    
#Pandoraea                                    
#Gottschalkia
#Desulfovibrio                                     
#Cyanothece                                 
#Macrococcus
#Desulfitobacterium
#Melittangium                                     
#Faecalibaculum
                                                 
g21 <- ggboxplot(G.tr.DF, x = "Status", y = "Denitrobacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Denitrobacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
g21 <- g21 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g21 <- g21 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g21 <- g21 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g22 <- ggboxplot(G.tr.DF, x = "Status", y = "Gemmatimonas", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Gemmatimonas"))) + theme(axis.text.x = element_blank()) + theme_classic()
g22 <- g22 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g22 <- g22 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g22 <- g22 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g23 <- ggboxplot(G.tr.DF, x = "Status", y = "Pandoraea", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Pandoraea"))) + theme(axis.text.x = element_blank()) + theme_classic()
g23 <- g23 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g23 <- g23 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g23 <- g23 + guides(fill=guide_legend(ncol=3,byrow=TRUE))
                                                 
g24 <- ggboxplot(G.tr.DF, x = "Status", y = "Gottschalkia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Gottschalkia"))) + theme(axis.text.x = element_blank()) + theme_classic()
g24 <- g24 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g24 <- g24 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g24 <- g24 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g25 <- ggboxplot(G.tr.DF, x = "Status", y = "Desulfovibrio", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Desulfovibrio"))) + theme(axis.text.x = element_blank()) + theme_classic()
g25 <- g25 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g25 <- g25 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g25 <- g25 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

                                 
g26 <- ggboxplot(G.tr.DF, x = "Status", y = "Cyanothece", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Cyanothece"))) + theme(axis.text.x = element_blank()) + theme_classic()
g26 <- g26 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g26 <- g26 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g26 <- g26 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g27 <- ggboxplot(G.tr.DF, x = "Status", y = "Macrococcus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Macrococcus"))) + theme(axis.text.x = element_blank()) + theme_classic()
g27 <- g27 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g27 <- g27 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g27 <- g27 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g28 <- ggboxplot(G.tr.DF, x = "Status", y = "Desulfitobacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Desulfitobacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
g28 <- g28 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g28 <- g28 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g28 <- g28 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g29 <- ggboxplot(G.tr.DF, x = "Status", y = "Melittangium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Melittangium"))) + theme(axis.text.x = element_blank()) + theme_classic()
g29 <- g29 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g29 <- g29 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g29 <- g29 + guides(fill=guide_legend(ncol=3,byrow=TRUE))

g30 <- ggboxplot(G.tr.DF, x = "Status", y = "Faecalibaculum", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Faecalibaculum"))) + theme(axis.text.x = element_blank()) + theme_classic()
g30 <- g30 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 5.5),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
g30 <- g30 + scale_fill_manual(labels = SL, values = g_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
g30 <- g30 + guides(fill=guide_legend(ncol=3,byrow=TRUE))
                                
tiff(filename = "KinshasaKonzo3_Genus_RF_ULPZ_Geo_Boxplots.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(s,s,g21, g22, g23, g24, g25, g26, g27, g28, g29, g30, labels = c("","","A","B","C","D","E","F","G","H","I","J"), font.label = list(size = 6), ncol = 4, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()  

#log 10 phyloseq                                                 
Control.G.tr.log10 <- transform_sample_counts(Control.G.tr, function(x) log10(x))
Disease.G.tr.log10 <- transform_sample_counts(Disease.G.tr, function(x) log10(x))     
                                                 
                                                 
#Supp 7: ULPZ vs. UHPZ
#ULPZ vs. UHPZ (Control)
                                                 
#Gordonibacter                                                 
#Denitrobacterium 
#Tumebacillus
#Adlercreutzia 
#Photobacterium
#Colwellia                                                 
#Slackia                                    
#Shewanella                                     
#Moraxella
#Tolumonas
                                     
my_comparisons <- list( c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone")) 

                                     
G <- Control.G.tr.log10
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Control.G.tr.log10@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

c1 <- ggboxplot(G.tr.DF, x = "Status", y = "Gordonibacter", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Gordonibacter"))) + theme(axis.text.x = element_blank()) + theme_classic()
c1 <- c1 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c1 <- c1 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c1 <- c1 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c2 <- ggboxplot(G.tr.DF, x = "Status", y = "Denitrobacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Denitrobacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
c2 <- c2 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c2 <- c2 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c2 <- c2 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c3 <- ggboxplot(G.tr.DF, x = "Status", y = "Tumebacillus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Tumebacillus"))) + theme(axis.text.x = element_blank()) + theme_classic()
c3 <- c3 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c3 <- c3 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c3 <- c3 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c4 <- ggboxplot(G.tr.DF, x = "Status", y = "Adlercreutzia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Adlercreutzia"))) + theme(axis.text.x = element_blank()) + theme_classic()
c4 <- c4 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c4 <- c4 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c4 <- c4 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c5 <- ggboxplot(G.tr.DF, x = "Status", y = "Photobacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Photobacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
c5 <- c5 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c5 <- c5 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c5 <- c5 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
c6 <- ggboxplot(G.tr.DF, x = "Status", y = "Colwellia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Colwellia"))) + theme(axis.text.x = element_blank()) + theme_classic()
c6 <- c6 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c6 <- c6 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c6 <- c6 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c7 <- ggboxplot(G.tr.DF, x = "Status", y = "Slackia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Slackia"))) + theme(axis.text.x = element_blank()) + theme_classic()
c7 <- c7 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c7 <- c7 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c7 <- c7 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c8 <- ggboxplot(G.tr.DF, x = "Status", y = "Shewanella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Shewanella"))) + theme(axis.text.x = element_blank()) + theme_classic()
c8 <- c8 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c8 <- c8 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c8 <- c8 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c9 <- ggboxplot(G.tr.DF, x = "Status", y = "Moraxella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Moraxella"))) + theme(axis.text.x = element_blank()) + theme_classic()
c9 <- c9 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c9 <- c9 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c9 <- c9 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

c10 <- ggboxplot(G.tr.DF, x = "Status", y = "Tolumonas", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Tolumonas"))) + theme(axis.text.x = element_blank()) + theme_classic()
c10 <- c10 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
c10 <- c10 + scale_fill_manual(labels = SSL, values = control_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2) #+ stat_compare_means( label.y = 1, size = 1)
c10 <- c10 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
tiff(filename = "KinshasaKonzo3_Genus_RF_ULPZ_UHPZ_Boxplots.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(s,s,c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, ncol = 4, nrow = 3, labels = c("","","A","B","C","D","E","F","G","H","I","J"), font.label = list(size = 6), common.legend = TRUE, legend = "bottom")
dev.off()    
                                                 
##Supp 9: KLPZ vs. KHPZ
#KLPZ vs. KHPZ (Disease)
                                     
#Adlercreutzia
#Herminiimonas
#Brevibacillus                                              
#Eggerthella                                  
#Flammeovirga 
#Escherichia
#Denitrobacterium                                                 
#Psychrobacter
#Acetoanaerobium                                                                
#Thermacetogenium                                  

                                     
my_comparisons <- list( c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone")) 

                                     
G <- Disease.G.tr.log10
                                               
G.tr_META <- as.data.frame(G@sam_data)
G.tr_OTU <- as.data.frame(t(G@otu_table))
G.tr.DF <- cbind(G.tr_OTU, G.tr_META$Status)

colnames(G.tr.DF)[colnames(G.tr.DF)=="G.tr_META$Status"] <- "Status"
for (i in 1:nrow(G.tr.DF))
  {G.tr.DF[i,]$Status <- Disease.G.tr.log10@sam_data[rownames(G.tr.DF[i,]),]$Status
  }
G.tr.DF$Status <- factor(G.tr.DF$Status, levels = c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
                                     
d1 <- ggboxplot(G.tr.DF, x = "Status", y = "Adlercreutzia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Adlercreutzia"))) + theme(axis.text.x = element_blank()) + theme_classic()
d1 <- d1 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d1 <- d1 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d1 <- d1 + guides(fill=guide_legend(ncol=2,byrow=TRUE))

d2 <- ggboxplot(G.tr.DF, x = "Status", y = "Herminiimonas", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Herminiimonas"))) + theme(axis.text.x = element_blank()) + theme_classic()
d2 <- d2 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d2 <- d2 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d2 <- d2 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
d3 <- ggboxplot(G.tr.DF, x = "Status", y = "Brevibacillus", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Brevibacillus"))) + theme(axis.text.x = element_blank()) + theme_classic()
d3 <- d3 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d3 <- d3 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d3 <- d3 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                                 
d4 <- ggboxplot(G.tr.DF, x = "Status", y = "Eggerthella", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Eggerthella"))) + theme(axis.text.x = element_blank()) + theme_classic()
d4 <- d4 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d4 <- d4 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d4 <- d4 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
d5 <- ggboxplot(G.tr.DF, x = "Status", y = "Flammeovirga", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Flammeovirga"))) + theme(axis.text.x = element_blank()) + theme_classic()
d5 <- d5 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d5 <- d5 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d5 <- d5 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                                 
d6 <- ggboxplot(G.tr.DF, x = "Status", y = "Escherichia", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Escherichia"))) + theme(axis.text.x = element_blank()) + theme_classic()
d6 <- d6 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d6 <- d6 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d6 <- d6 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
d7 <- ggboxplot(G.tr.DF, x = "Status", y = "Denitrobacterium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Denitrobacterium"))) + theme(axis.text.x = element_blank()) + theme_classic()
d7 <- d7 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d7 <- d7 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d7 <- d7 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
d8 <- ggboxplot(G.tr.DF, x = "Status", y = "Psychrobacter", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Psychrobacter"))) + theme(axis.text.x = element_blank()) + theme_classic()
d8 <- d8 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d8 <- d8 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d8 <- d8 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                     
d9 <- ggboxplot(G.tr.DF, x = "Status", y = "Acetoanaerobium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Acetoanaerobium"))) + theme(axis.text.x = element_blank()) + theme_classic()
d9 <- d9 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d9 <- d9 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d9 <- d9 + guides(fill=guide_legend(ncol=2,byrow=TRUE))
                                                                         
d10 <- ggboxplot(G.tr.DF, x = "Status", y = "Thermacetogenium", fill = "Status", xlab = "Samples", ylab = "log(rel. abund.)", title = expression(italic("Thermacetogenium"))) + theme(axis.text.x = element_blank()) + theme_classic()
d10 <- d10 + theme(legend.position="right") + theme(legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8), legend.title = element_blank(), legend.border = NULL) + guides(fill=guide_legend(ncol=1,byrow=TRUE)) + theme(title = element_text(size = 6),axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), axis.text.y = element_text(size = 8)) + stat_boxplot(geom ='errorbar')
d10 <- d10 + scale_fill_manual(labels = SSL, values = disease_color) + stat_compare_means(comparisons = my_comparisons, label = "p.signif", method = "wilcox.test", size = 2)
d10 <- d10 + guides(fill=guide_legend(ncol=2,byrow=TRUE))                                   

tiff(filename = "KinshasaKonzo3_Genus_RF_KLPZ_KHPZ_Boxplots.tiff", width = 7, height = 7, units = "in", res = 600)
ggarrange(s,s, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, labels = c("","","A","B","C","D","E","F","G","H","I","J"), font.label = list(size = 6), ncol = 4, nrow = 3, common.legend = TRUE, legend = "bottom")
dev.off()   
                                              
                                              
                                              
                                             
#------------------------FUNTIONAL---------------------------------
                                              
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostBracken/KinshasaControl_Konzo3_PostBracken")

#META
Konzo_meta <- read.csv("./KinshasaControl_Konzo3_Meta_Mod.csv")
names(Konzo_meta)<-c("Sample","Name","Run","ID","Region","Status","Disease","Sample_ID","Collection_date","DNA_Concentration","Isolation_date","Elution","Age","Sex","Disease_Old","Geography")
rownames(Konzo_meta)<-as.character(Konzo_meta[,1])
META<-sample_data(Konzo_meta)


#Konzo Functional KO 
setwd("~/Dropbox/Konzo_Microbiome/Konzo1Konzo3/Konzo1_Konzo3_PostAlbanFunctional")

#KO count
Konzo_ko_count <- read.csv("./KinshasaKonzo3_KO_Count.csv")
Konzo_ko <- read.csv("./KinshasaKonzo3_KO.csv")

Konzo_KO_count <-as.matrix(unname(Konzo_ko_count[1:nrow(Konzo_ko_count),2:(ncol(Konzo_ko_count))]))
rownames(Konzo_KO_count)<-as.character( Konzo_ko_count[,1])
nam <-names(Konzo_ko_count)
colnames(Konzo_KO_count)<-c(as.character(nam[2:length(nam)]))
KO_count = otu_table(Konzo_KO_count, taxa_are_rows = TRUE)

#TAX (KO)
Konzo_KO<-as.matrix(unname(Konzo_ko[,2]))
rownames(Konzo_KO)<-as.character(unname(Konzo_KO[,1]))
colnames(Konzo_KO)<-"KO"
KO = tax_table(Konzo_KO)

#Create the PhyloseqObject
KonzoData_KO_count <-phyloseq(KO_count, KO, META)

#Set NAs to 0
KonzoData_KO_count@otu_table[is.na(KonzoData_KO_count@otu_table)] <- 0
KonzoData_KO_count@sam_data$Status <- factor(KonzoData_KO_count@sam_data$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))
#write.csv(KonzoData_KO_count@otu_table, file = "./KinshasaKonzo3_KO_Count_BySample.csv")

#Merge samples by group/status                                         
KonzoData_KO_count_status <- merge_samples(KonzoData_KO_count, KonzoData_KO_count@sam_data$Status) #merge_smaples by default sums the values for otu
KonzoData_KO_count_status <- transform_sample_counts(KonzoData_KO_count_status, function(x) x / 30) #average the sum of relabund in each group                                                                                                                                    
#write.csv(t(KonzoData_KO_count_status@otu_table), file = "./KinshasaKonzo3_KO_AvgCount_ByStatus.csv")
                                                  
#Read Counts to Relative Abundance
KonzoData_KO_tr <- transform_sample_counts(KonzoData_KO_count, function(x) x / sum(x))                                  
#write.csv(KonzoData_KO_tr@otu_table, file = "./KinshasaKonzo3_KO_RelAbund.csv")  

#Merge samples by group/status                                         
KonzoData_KO_tr_status <- merge_samples(KonzoData_KO_tr, KonzoData_KO_tr@sam_data$Status) #merge_smaples by default sums the values for otu
KonzoData_KO_tr_status <- transform_sample_counts(KonzoData_KO_tr_status, function(x) x / 30) #average the sum of relabund in each group                                                                                                                                    
#write.csv(t(KonzoData_KO_tr_status@otu_table), file = "./KinshasaKonzo3_KO_AvgRelAbund_ByStatus.csv")
                                                 
# Filtering
Kinshasa.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Kinshasa", KonzoData_KO_tr)
Masimanimba.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Masimanimba", KonzoData_KO_tr)
ULPZ.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Unaffected_Low_Prevalence_Zone", KonzoData_KO_tr)
KLPZ.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData_KO_tr)
UHPZ.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData_KO_tr)
KHPZ.KO.tr <- prune_samples(KonzoData_KO_tr@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData_KO_tr)
                                     
Kinshasa.KO.tr.f <- filter_taxa(Kinshasa.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
Masimanimba.KO.tr.f <- filter_taxa(Masimanimba.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
ULPZ.KO.tr.f <- filter_taxa(ULPZ.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KLPZ.KO.tr.f <- filter_taxa(KLPZ.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
UHPZ.KO.tr.f <- filter_taxa(UHPZ.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)
KHPZ.KO.tr.f <- filter_taxa(KHPZ.KO.tr, function (x) mean(x) >= 1e-4, prune = TRUE)

filterList1 <- union(Kinshasa.KO.tr.f@tax_table,Masimanimba.KO.tr.f@tax_table) #Kin, Mas
filterList2 <- union(ULPZ.KO.tr.f@tax_table, KLPZ.KO.tr.f@tax_table) #ULPZ, KLPZ
filterList3 <- union(UHPZ.KO.tr.f@tax_table,KHPZ.KO.tr.f@tax_table)
filterList4 <- union(filterList1, filterList2) #Kin, Mas, ULPZ, KLPZ
filterList <- union(filterList3,filterList4) # Kin, Mas, ULPS, KLPZ,UHPZ, KHPZ

#write.csv(filterList, file = "Kinshasa_Konzo3_KO_f_0.0001.csv")
                            
x <- read.csv("Kinshasa_Konzo3_KO_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                                                 
KonzoData_KO_count.f <- prune_taxa(f_0.0001, KonzoData_KO_count) #filtered readcount phyloseq object
KonzoData_KO_tr.f <- prune_taxa(f_0.0001, KonzoData_KO_tr) #filtered rel abund phyloseq object                                            
KonzoData_KO_tr_status.f <- prune_taxa(f_0.0001,KonzoData_KO_tr_status)
                            
KonzoData_KO_tr.f.otu <- as.data.frame(KonzoData_KO_tr.f@otu_table)                           
KonzoData_KO_tr_status.f.otu <- as.data.frame(t(KonzoData_KO_tr_status.f@otu_table))                            
 
KO_func_f_0.0001 <- read.csv("./KinshasaKonzo3_KO_f_0.0001_function_AM.csv", row.names = 1) 

                                                      
KonzoData_KO_tr.f.otu.func <- cbind(KO_func_f_0.0001, KonzoData_KO_tr.f.otu)
#write.csv(KonzoData_KO_tr.f.otu.func, file = "./KinshasaKonzo3_KO_f_0.0001_RelAbund_WithFunction.csv")
                            
KonzoData_KO_tr_status.f.otu.func <- cbind(KO_func_f_0.0001, KonzoData_KO_tr_status.f.otu)
#write.csv(KonzoData_KO_tr_status.f.otu.func, file = "./KinshasaKonzo3_KO_f_0.0001_AvgRelAbund_ByStatus_WithFunction.csv")                                             

                            
#Mean and Standard Deviation
KonzoData_KO_tr.df <- as.data.frame(t(KonzoData_KO_tr@otu_table))
KonzoData_KO_tr.df <- cbind(KonzoData_KO_tr.df, KonzoData_KO_tr@sam_data$Status)

colnames(KonzoData_KO_tr.df)[colnames(KonzoData_KO_tr.df)=="KonzoData_KO_tr@sam_data$Status"] <- "Status"
for (i in 1:nrow(KonzoData_KO_tr.df))
  {KonzoData_KO_tr.df[i,]$Status <- KonzoData_KO_tr@sam_data[rownames(KonzoData_KO_tr.df[i,]),]$Status
  } 
                                                 
KonzoData_KO_tr.avg <- KonzoData_KO_tr.df %>% group_by(Status) %>% summarise_each(funs(mean)) 
KonzoData_KO_tr.avg.x <- t(KonzoData_KO_tr.avg)
colnames(KonzoData_KO_tr.avg.x) <- KonzoData_KO_tr.avg.x[1,]   
KonzoData_KO_tr.avg.x <- KonzoData_KO_tr.avg.x[-1,]
colnames(KonzoData_KO_tr.avg.x) <- paste("Avg", colnames(KonzoData_KO_tr.avg.x), sep = "_")                                                 
                                                
                                                 
KonzoData_KO_tr.sd <- KonzoData_KO_tr.df %>% group_by(Status) %>% summarise_each(funs(sd))  
KonzoData_KO_tr.sd.x <- t(KonzoData_KO_tr.sd)
colnames(KonzoData_KO_tr.sd.x) <- KonzoData_KO_tr.sd.x[1,]   
KonzoData_KO_tr.sd.x <- KonzoData_KO_tr.sd.x[-1,]
colnames(KonzoData_KO_tr.sd.x) <- paste("SD", colnames(KonzoData_KO_tr.sd.x), sep = "_")                                                 

KonzoData_KO_tr.avg.sd <-merge(KonzoData_KO_tr.avg.x, KonzoData_KO_tr.sd.x,by='row.names', sort = FALSE) 
rownames(KonzoData_KO_tr.avg.sd) <- KonzoData_KO_tr.avg.sd[,1]   
KonzoData_KO_tr.avg.sd <- KonzoData_KO_tr.avg.sd[,-1]  
                           
KonzoData_KO_tr.avg.sd <- KonzoData_KO_tr.avg.sd[, c(1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12)]                                                   
write.csv(KonzoData_KO_tr.avg.sd, file = "./KonzoData_KO_AvgRelAbund_SD_ByGroup.csv") 
                           
KonzoData_KO_tr.avg.sd.f <- subset(KonzoData_KO_tr.avg.sd, rownames(KonzoData_KO_tr.avg.sd) %in% f_0.0001)                                             
write.csv(KonzoData_KO_tr.avg.sd.f, file = "./KonzoData_KO_AvgRelAbund_SD_ByGroup_filtered.csv") 
                                                     
KO.tr <- merge(KonzoData_KO_tr.avg.sd,as.data.frame(KonzoData_KO_tr@otu_table),by='row.names', sort = FALSE)                           
write.csv(KO.tr, file = "./KonzoData_KO_RelAbund_Supp.csv")      
                            
KO_func_f_0.0001 <- read.csv("./KinshasaKonzo3_KO_f_0.0001_function_AM.csv", row.names = 1) 
                                                                                
KO.tr.f <- merge(KonzoData_KO_tr.avg.sd.f,as.data.frame(KonzoData_KO_tr.f@otu_table),by='row.names', sort = FALSE)  
rownames(KO.tr.f) <- KO.tr.f[,1]   
KO.tr.f <- KO.tr.f[,-1]                              
                            
KO.tr.f_func <- cbind(KO_func_f_0.0001, KO.tr.f)
                            
write.csv(KO.tr.f_func, file = "./KonzoData_KO_RelAbund_Func_filtered_Supp.csv")                            

                            
#Initially perform a multivariate analysis to determine which variables influence the gut microbiome of the various cohorts    
                            
KonzoData_KO_tr.sam <- as.data.frame(as.matrix(sample_data(KonzoData_KO_tr)))
KonzoData_KO_tr.sam$Status <- as.factor(KonzoData_KO_tr.sam$Status)
KonzoData_KO_tr.sam$Status <- factor(KonzoData_KO_tr.sam$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(KonzoData_KO_tr, method="bray")
bdiv_bray <- adonis(brayd ~ KonzoData_KO_tr.sam$Geography * KonzoData_KO_tr.sam$Region * KonzoData_KO_tr.sam$Disease * KonzoData_KO_tr.sam$Age * KonzoData_KO_tr.sam$Sex, perm=99999); bdiv_bray
                            
                            
#Supplementary Figure 3                            
#Geography_KO  
                            
x <- read.csv("Kinshasa_Konzo3_KO_f_0.0001.csv", row.names = 1, colClasses = "character")
f_0.0001 <- unlist(x)
                            
Geography.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status != "Konzo_Low_Prevalence_Zone") & (KonzoData_KO_tr@sam_data$Status != "Konzo_High_Prevalence_Zone"), KonzoData_KO_tr)                                              

Geography.KO.tr.f <- prune_taxa(f_0.0001, Geography.KO.tr) 
                                                                                                                                                                                                                                                                
Geography.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(Geography.KO.tr.f)))
Geography.KO.tr.f.sam$Status <- as.factor(Geography.KO.tr.f.sam$Status)                       
Geography.KO.tr.f.sam$Status <- factor(Geography.KO.tr.f.sam$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

                            
brayd <- phyloseq::distance(Geography.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ Geography.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_Geography_KO.tr.filtered.txt") #1 e-5
                                                  
                                                  
p1 = plot_ordination(Geography.KO.tr.f, ordinate(Geography.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 2, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
ko_PGB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = geography_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA), legend.margin=margin(c(-8,0,0,0))) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))


ko_PGBt <- ko_PGB + stat_ellipse(type = "t") + scale_x_continuous(position = "top") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
ko_PGBt <- ko_PGBt + theme(legend.position="bottom")
ko_PGBt <- ko_PGBt + annotate("text", x = -0.24, y = -0.39, label = expression(paste("p = 1x",10^-5)), size = 2.5)
ko_PGBt <- ggarrange(ko_PGBt,labels = c("A"),font.label = list(size = 7))
                                                 

#KinMas                                                   
KinMas.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Kinshasa" | KonzoData_KO_tr@sam_data$Status == "Masimanimba"), KonzoData_KO_tr)  
                            
KinMas.KO.tr.f <- prune_taxa(f_0.0001, KinMas.KO.tr) 
                            
KinMas.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(KinMas.KO.tr.f)))
KinMas.KO.tr.f.sam$Status <- as.factor(KinMas.KO.tr.f.sam$Status)
KinMas.KO.tr.f.sam$Status <- factor(KinMas.KO.tr.f.sam$Status, levels = c("Kinshasa", "Masimanimba"))

brayd <- phyloseq::distance(KinMas.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ KinMas.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_KinMas_KO.tr.filtered.txt") #0.00093                                                 
                                                  
p1 = plot_ordination(KinMas.KO.tr.f, ordinate(KinMas.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
                                                  
ko_PKMB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinmas_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

ko_PKMBt <- ko_PKMB + stat_ellipse(type = "t") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
ko_PKMBt <- ko_PKMBt + theme(legend.position="none")
ko_PKMBt <- ko_PKMBt + annotate("text", x = -0.31, y = -0.3, label = expression(paste("p = 0.00093")), size = 2)
ko_PKMBt <- ggarrange(ko_PKMBt,labels = c("B"),font.label = list(size = 7))                                                  
                                                  

#KinULPZ
KinULPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Kinshasa" | KonzoData_KO_tr@sam_data$Status == "Unaffected_Low_Prevalence_Zone"), KonzoData_KO_tr)  
                            
KinULPZ.KO.tr.f <- prune_taxa(f_0.0001, KinULPZ.KO.tr) 

                            
KinULPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(KinULPZ.KO.tr.f)))
KinULPZ.KO.tr.f.sam$Status <- as.factor(KinULPZ.KO.tr.f.sam$Status)
KinULPZ.KO.tr.f.sam$Status <- factor(KinULPZ.KO.tr.f.sam$Status, levels = c("Kinshasa", "Unaffected_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(KinULPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ KinULPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_KinULPZ_KO.tr.filtered.txt")      #0.00058                                            

p1 = plot_ordination(KinULPZ.KO.tr.f, ordinate(KinULPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]

ko_PKUB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinulpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

ko_PKUBt <- ko_PKUB + stat_ellipse(type = "t") + theme(plot.margin=unit(c(0.15,0.15,0.15,0.15), "lines"))
ko_PKUBt <- ko_PKUBt + theme(legend.position="none")
ko_PKUBt <- ko_PKUBt + annotate("text", x = -0.17, y = -0.22, label = expression(paste("p = 0.00058")), size = 2)
ko_PKUBt <- ggarrange(ko_PKUBt,labels = c("C"),font.label = list(size = 7))                                                  
                                                  
                                                  
#KinUHPZ  
KinUHPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Kinshasa" | KonzoData_KO_tr@sam_data$Status == "Unaffected_High_Prevalence_Zone"), KonzoData_KO_tr) 
KinUHPZ.KO.tr.f <- prune_taxa(f_0.0001, KinUHPZ.KO.tr)                             
                            
KinUHPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(KinUHPZ.KO.tr.f)))
KinUHPZ.KO.tr.f.sam$Status <- as.factor(KinUHPZ.KO.tr.f.sam$Status)
KinUHPZ.KO.tr.f.sam$Status <- factor(KinUHPZ.KO.tr.f.sam$Status, levels = c("Kinshasa", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(KinUHPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ KinUHPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_KinUHPZ_KO.tr.filtered.txt")     #2 e -5                                             

p1 = plot_ordination(KinUHPZ.KO.tr.f, ordinate(KinUHPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
                                                
ko_PKUHB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = kinuhpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

ko_PKUHBt <- ko_PKUHB + stat_ellipse(type = "t") +  scale_y_continuous(position = "right")+ theme(plot.margin=unit(c(0.15,0.15,0.15,0.6), "lines"))
ko_PKUHBt <- ko_PKUHBt + theme(legend.position="none")                                                                                                   
ko_PKUHBt <- ko_PKUHBt + annotate("text", x = -0.17, y = -0.22, label = expression(paste("p = 2x",10^-5)), size = 2)
ko_PKUHBt <- ggarrange(ko_PKUHBt,labels = c("D"),font.label = list(size = 7))
                                                 
#MasULPZ
MasULPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Masimanimba" | KonzoData_KO_tr@sam_data$Status == "Unaffected_Low_Prevalence_Zone"), KonzoData_KO_tr)   
MasULPZ.KO.tr.f <- prune_taxa(f_0.0001, MasULPZ.KO.tr)  
                            
MasULPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(MasULPZ.KO.tr.f)))
MasULPZ.KO.tr.f.sam$Status <- as.factor(MasULPZ.KO.tr.f.sam$Status)
MasULPZ.KO.tr.f.sam$Status <- factor(MasULPZ.KO.tr.f.sam$Status, levels = c("Masimanimba", "Unaffected_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(MasULPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ MasULPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_MasULPZ_KO.tr.filtered.txt")       #5e-05                                          

p1 = plot_ordination(MasULPZ.KO.tr.f, ordinate(MasULPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
ko_PMUB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = masulpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

ko_PMUBt <- ko_PMUB + stat_ellipse(type = "t") + scale_x_continuous(position = "top") + scale_y_continuous(position = "right") + theme(plot.margin=unit(c(0.15,0.15,0.25,0.25), "lines"))
ko_PMUBt <- ko_PMUBt + theme(legend.position="none")
ko_PMUBt <- ko_PMUBt + annotate("text", x = -0.21, y = -0.27, label = expression(paste("p = 5x",10^-5)), size = 2)
ko_PMUBt <- ggarrange(ko_PMUBt,labels = c("E"),font.label = list(size = 7))
                                                  
#MasUHPZ  
MasUHPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Masimanimba" | KonzoData_KO_tr@sam_data$Status == "Unaffected_High_Prevalence_Zone"), KonzoData_KO_tr)  
                            
MasUHPZ.KO.tr.f <- prune_taxa(f_0.0001, MasUHPZ.KO.tr)  
                            
MasUHPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(MasUHPZ.KO.tr.f)))
MasUHPZ.KO.tr.f.sam$Status <- as.factor(MasUHPZ.KO.tr.f.sam$Status)
MasUHPZ.KO.tr.f.sam$Status <- factor(MasUHPZ.KO.tr.f.sam$Status, levels = c("Masimanimba", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(MasUHPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ MasUHPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_MasUHPZ_KO.tr.filtered.txt")   #1e-04                                               

p1 = plot_ordination(MasUHPZ.KO.tr.f, ordinate(MasUHPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
#p1 <- as_ggplot(p1)
ko_PMUHB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = masuhpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=7)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7))

ko_PMUHBt <- ko_PMUHB + stat_ellipse(type = "t") + scale_x_continuous(position = "top") + scale_y_continuous(position = "right") + theme(plot.margin=unit(c(0.15,0.15,0.85,0.25), "lines"))
ko_PMUHBt <- ko_PMUHBt + theme(legend.position="none")
ko_PMUHBt <- ko_PMUHBt + annotate("text", x = -0.22, y = -0.34, label = expression(paste("p = 1x", 10^-4)), size = 2)
ko_PMUHBt <- ggarrange(ko_PMUHBt,labels = c("F"),font.label = list(size = 7))
                                                  
Geo <- arrangeGrob(ko_PGBt, ko_PMUBt, ko_PMUHBt, ko_PKMBt, ko_PKUBt, ko_PKUHBt,                             
             ncol = 6, nrow = 6,
             layout_matrix = rbind(c(1,1,1,1,2,2), c(1,1,1,1,2,2), c(1,1,1,1,3,3), c(1,1,1,1,3,3), c(4,4,5,5,6,6), c(4,4,5,5,6,6)))

tiff(filename = "Geography_AllUnaffected_KO_Filtered_PCoA.tiff", width = 5.5, height = 5.5, units = "in", res = 600)
ggarrange(as_ggplot(Geo))
dev.off()

                            
                            
#KinKLPZ
KinKLPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Kinshasa" | KonzoData_KO_tr@sam_data$Status == "Konzo_Low_Prevalence_Zone"), KonzoData_KO_tr)                                                                                                                                                                                                                                                                                                             
KinKLPZ.KO.tr.f <- prune_taxa(f_0.0001, KinKLPZ.KO.tr)
                            
KinKLPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(KinKLPZ.KO.tr.f)))
KinKLPZ.KO.tr.f.sam$Status <- as.factor(KinKLPZ.KO.tr.f.sam$Status)
KinKLPZ.KO.tr.f.sam$Status <- factor(KinKLPZ.KO.tr.f.sam$Status, levels = c("Kinshasa", "Konzo_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(KinKLPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ KinKLPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_KinKLPZ_KO.tr.filtered.txt") #0.0031                                                 

#KinKHPZ
KinKHPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Kinshasa" | KonzoData_KO_tr@sam_data$Status == "Konzo_High_Prevalence_Zone"), KonzoData_KO_tr)  
KinKHPZ.KO.tr.f <- prune_taxa(f_0.0001, KinKHPZ.KO.tr)                            
                            
                            
KinKHPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(KinKHPZ.KO.tr.f)))
KinKHPZ.KO.tr.f.sam$Status <- as.factor(KinKHPZ.KO.tr.f.sam$Status)
KinKHPZ.KO.tr.f.sam$Status <- factor(KinKHPZ.KO.tr.f.sam$Status, levels = c("Kinshasa", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(KinKHPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ KinKHPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_genus_adonis_KinKHPZ_KO.tr.filtered.txt")    #1e-05                                              
                                                  
#MasKLPZ
MasKLPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Masimanimba" | KonzoData_KO_tr@sam_data$Status == "Konzo_Low_Prevalence_Zone"), KonzoData_KO_tr)  
MasKLPZ.KO.tr.f <- prune_taxa(f_0.0001, MasKLPZ.KO.tr)                            
                            
                            
MasKLPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(MasKLPZ.KO.tr.f)))
MasKLPZ.KO.tr.f.sam$Status <- as.factor(MasKLPZ.KO.tr.f.sam$Status)
MasKLPZ.KO.tr.f.sam$Status <- factor(MasKLPZ.KO.tr.f.sam$Status, levels = c("Masimanimba", "Konzo_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(MasKLPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ MasKLPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_MasKLPZ_KO.tr.filtered.txt")    #0.00045                                             

#MasKHPZ
MasKHPZ.KO.tr <- prune_samples((KonzoData_KO_tr@sam_data$Status == "Masimanimba" | KonzoData_KO_tr@sam_data$Status == "Konzo_High_Prevalence_Zone"), KonzoData_KO_tr) 
MasKHPZ.KO.tr.f <- prune_taxa(f_0.0001, MasKHPZ.KO.tr)                            
                            
                            
MasKHPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(MasKHPZ.KO.tr.f)))
MasKHPZ.KO.tr.f.sam$Status <- as.factor(MasKHPZ.KO.tr.f.sam$Status)
MasKHPZ.KO.tr.f.sam$Status <- factor(MasKHPZ.KO.tr.f.sam$Status, levels = c("Masimanimba", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(MasKHPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ MasKHPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_MasKHPZ_KO.tr.filtered.txt")   #3e-05                                               
     
                            
#Supplementary Figure 8
                            
#Control                                                  
Control.KO.tr <-  prune_samples(KonzoData_KO_tr@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData_KO_tr@sam_data$Status == "Unaffected_High_Prevalence_Zone", KonzoData_KO_tr)
Control.KO.tr.f <- prune_taxa(f_0.0001, Control.KO.tr)                            
                                                       
                            
Control.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(Control.KO.tr.f)))
Control.KO.tr.f.sam$Status <- as.factor(Control.KO.tr.f.sam$Status)
Control.KO.tr.f.sam$Status <- factor(Control.KO.tr.f.sam$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"))

brayd <- phyloseq::distance(Control.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ Control.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_Control_KO.tr.filtered.txt")  #0.05741       

p1 = plot_ordination(Control.KO.tr.f, ordinate(Control.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
ko_PCB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = control_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

ko_PCBt <- ko_PCB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
ko_PCBt <- ko_PCBt + annotate("text", x = 0.25, y = -0.3, label = expression(paste("p = 0.05741")), size = 2)
                            
                            
                                            
#Disease                                                  
Disease.KO.tr <-  prune_samples(KonzoData_KO_tr@sam_data$Status == "Konzo_Low_Prevalence_Zone" | KonzoData_KO_tr@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData_KO_tr)
Disease.KO.tr.f <- prune_taxa(f_0.0001, Disease.KO.tr)                            
                                                        
Disease.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(Disease.KO.tr.f)))
Disease.KO.tr.f.sam$Status <- as.factor(Disease.KO.tr.f.sam$Status)
Disease.KO.tr.f.sam$Status <- factor(Disease.KO.tr.f.sam$Status, levels = c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(Disease.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ Disease.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_Disease_KO.tr.filtered.txt")     #0.0514                                          

p1 = plot_ordination(Disease.KO.tr.f, ordinate(Disease.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
ko_PKB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = disease_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

ko_PKBt <- ko_PKB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
ko_PKBt <- ko_PKBt + annotate("text", x = 0.22, y = -0.3, label = expression(paste("p = 0.0514")), size = 2)
     
                            
#LPZ                                                  
LPZ.KO.tr <-  prune_samples(KonzoData_KO_tr@sam_data$Status == "Unaffected_Low_Prevalence_Zone" | KonzoData_KO_tr@sam_data$Status == "Konzo_Low_Prevalence_Zone", KonzoData_KO_tr)
LPZ.KO.tr.f <- prune_taxa(f_0.0001, LPZ.KO.tr)                            
                                                        
                            
LPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(LPZ.KO.tr.f)))
LPZ.KO.tr.f.sam$Status <- as.factor(LPZ.KO.tr.f.sam$Status)
LPZ.KO.tr.f.sam$Status <- factor(LPZ.KO.tr.f.sam$Status, levels = c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"))

brayd <- phyloseq::distance(LPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ LPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_LPZ_KO.tr.filtered.txt")    #0.8929                                          

p1 = plot_ordination(LPZ.KO.tr.f, ordinate(LPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
ko_PNIB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = lpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

ko_PNIBt <- ko_PNIB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
ko_PNIBt <- ko_PNIBt + annotate("text", x = 0.13, y = -0.3, label = expression(paste("p = 0.8929")), size = 2)
                                                  
#HPZ                                                 
HPZ.KO.tr <-  prune_samples(KonzoData_KO_tr@sam_data$Status == "Unaffected_High_Prevalence_Zone" | KonzoData_KO_tr@sam_data$Status == "Konzo_High_Prevalence_Zone", KonzoData_KO_tr)
HPZ.KO.tr.f <- prune_taxa(f_0.0001, HPZ.KO.tr)                            
                            
                            
HPZ.KO.tr.f.sam <- as.data.frame(as.matrix(sample_data(HPZ.KO.tr.f)))
HPZ.KO.tr.f.sam$Status <- as.factor(HPZ.KO.tr.f.sam$Status)
HPZ.KO.tr.f.sam$Status <- factor(HPZ.KO.tr.f.sam$Status, levels = c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

brayd <- phyloseq::distance(HPZ.KO.tr.f, method="bray")
bdiv_bray <- adonis(brayd ~ HPZ.KO.tr.f.sam$Status, perm=99999); bdiv_bray
#capture.output(bdiv_bray, file="relabund_bdiv_adonis_HPZ_KO.tr.filtered.txt")  #0.8634                                            

p1 = plot_ordination(HPZ.KO.tr.f, ordinate(HPZ.KO.tr.f, method="PCoA", dist="bray"), type="samples", color="Status") +
  geom_point(size = 1, stroke = 0, shape = 16)
p1$layers <- p1$layers[-1]
ko_PIB <- p1 + 
  labs(color = "Groups")+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
                                axis.line = element_line(colour = "black")) + scale_color_manual(values = hpz_color, labels = SSSL)+
  theme(legend.title=element_blank(), legend.margin=margin(-5,0,0,0), legend.position = "bottom", legend.background = element_rect(colour = NA, fill = NA)) + theme (legend.key = element_rect(colour = NA, fill = NA ), panel.border = element_rect(colour = "black", fill=NA, size=0.5)) + theme(legend.key.size = unit(.1, "cm")) + theme(legend.text = element_text(size=5)) +
  theme(axis.title.y = element_text(size = 7), axis.title.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6))

ko_PIBt <- ko_PIB + stat_ellipse(type = "t") + guides(fill=guide_legend(nrow=1))
ko_PIBt <- ko_PIBt + annotate("text", x = 0.21, y = -0.3, label = expression(paste("p = 0.8634")), size = 2)
                            
tiff(filename = "Kahemba_KO_Filtered_Control_Disease_LPZ_HPZ_PCoA.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
ggarrange(ko_PCBt, ko_PKBt, ko_PNIBt, ko_PIBt, labels = c("A","B", "C", "D"), ncol = 2, nrow = 2, font.label = list(size = 7))
dev.off()                                                  
                                                  
           
#Figure 6:LAB, Beta-glucosidase and Rhodanese

S <- KonzoData.S.tr

S.tr_META <- as.data.frame(S@sam_data)
S.tr_OTU <- as.data.frame(t(S@otu_table))
S.tr.DF <- cbind(S.tr_OTU, S.tr_META$Status)
S.tr.DF <- cbind(S.tr.DF, S.tr_META$Geography)

colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Status"] <- "Status"
colnames(S.tr.DF)[colnames(S.tr.DF)=="S.tr_META$Geography"] <- "Geography"

for (i in nrow(S.tr.DF))
{S.tr.DF[i,]$Status <- KonzoData.S.tr@sam_data[rownames(S.tr.DF[i,]),]$Status
}
S.tr.DF$Status <- factor(S.tr.DF$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

for (i in nrow(S.tr.DF))
{S.tr.DF[i,]$Geography <- KonzoData.S.tr@sam_data[rownames(S.tr.DF[i,]),]$Geography
}
S.tr.DF$Geography <- factor(S.tr.DF$Geography, levels = c("Kinshasa", "Masimanimba", "Low_Prevalence_Zone", "High_Prevalence_Zone"))

#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone), c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone")) 


colnames(S.tr.DF)<- gsub( " ", ".", colnames(S.tr.DF)) 

S.tr.DF.status <- melt(S.tr.DF[,c('Status', 'Leuconostoc.mesenteroides', 'Lactococcus.lactis', 'Lactobacillus.plantarum')],id.vars = 1)
S.tr.DF.geography <- melt(S.tr.DF[,c('Geography','Leuconostoc.mesenteroides', 'Lactococcus.lactis', 'Lactobacillus.plantarum')],id.vars = 1)

# for box plot: facet_zoom(ylim = c(0, 0.0015))


#temp <- c(`Leuconostoc.mesenteroides` = "Leuconostoc mesenteroides", `Lactobacillus.plantarum` = "Lactobacillus plantarum", `Lactococcus.lactis` = "Lactococcus lactis")
temp <- c(`Leuconostoc.mesenteroides` = "L. mesenteroides", `Lactobacillus.plantarum` = "L. plantarum", `Lactococcus.lactis` = "L. lactis")

errors = aggregate(. ~ Status + variable, 
                   data = S.tr.DF.status, 
                   FUN = sd)
errors2 = aggregate(. ~ Status + variable, 
                   data = S.tr.DF.status, 
                   FUN = se)
means = aggregate(. ~ Status + variable, 
                  data = S.tr.DF.status, FUN = mean)
                            
t6 <- ggplot(S.tr.DF.status,aes(x = Status,y = value)) + 
    geom_boxplot(aes(fill = variable), lwd=0.2, outlier.size = 0.15, fatten = 0.6) + theme_classic() + ylab("rel. abund.")
t6 <- t6 + theme(legend.position="bottom", legend.margin=margin(0,0,0,0)) + scale_x_discrete(labels= SSSL) + theme(plot.title = element_blank(), legend.key.size = unit(.3, "cm"), legend.text = element_text(size = 6, face = "italic"), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7), axis.title.x = element_blank())
t6 <- t6 + scale_fill_discrete(labels = temp)
t6 <- t6 + coord_cartesian(ylim = c(0, 0.0026)) + scale_y_continuous(breaks= seq(0.0005, 0.0026, by = 0.001), expand = c(0,0))

                                    
t7 <- ggplot(S.tr.DF.status,aes(x = Status,y = value)) + 
    geom_boxplot(aes(fill = variable), lwd=0.2,outlier.size = 0.15, fatten = 0.6) + theme_classic() + ylab("") #  scale_x_discrete(labels= SSSL, position = "top")
t7 <- t7 + theme(legend.position="bottom", legend.margin=margin(0,0,0,0)) + scale_x_discrete(labels= SSSL, position = "top") + theme(plot.title = element_blank(), legend.key.size = unit(.3, "cm"), legend.text = element_text(size = 6, face = "italic"), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 7), axis.title.x = element_blank())
t7 <- t7 + scale_fill_discrete(labels = temp) 
t7 <- t7 + coord_cartesian(ylim = c(0.0035, 0.037)) + scale_y_continuous(breaks= seq(0.0035, 0.037, by = 0.01), expand = c(0,0))  
t7 <- t7 + theme(legend.position = "NA")
                                    
lab <- ggarrange(t7, t6, ncol = 1 ,heights = c(1,2), align = "v")                                
 
tiff(filename = "Kinshasa_Konzo3_LAB_Species_BoxPlot_split.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
lab                                    
dev.off()

###
#Beta-Glucosidase
                                    
K <- KonzoData_KO_tr

K.tr_META <- as.data.frame(K@sam_data)
K.tr_OTU <- as.data.frame(t(K@otu_table))
K.tr.DF <- cbind(K.tr_OTU, K.tr_META$Status)
K.tr.DF <- cbind(K.tr.DF, K.tr_META$Geography)

colnames(K.tr.DF)[colnames(K.tr.DF)=="K.tr_META$Status"] <- "Status"
colnames(K.tr.DF)[colnames(K.tr.DF)=="K.tr_META$Geography"] <- "Geography"

for (i in nrow(K.tr.DF))
{K.tr.DF[i,]$Status <- KonzoData_KO_tr@sam_data[rownames(K.tr.DF[i,]),]$Status
}
K.tr.DF$Status <- factor(K.tr.DF$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

for (i in nrow(K.tr.DF))
{K.tr.DF[i,]$Geography <- KonzoData_KO_tr@sam_data[rownames(K.tr.DF[i,]),]$Geography
}
K.tr.DF$Geography <- factor(K.tr.DF$Geography, levels = c("Kinshasa", "Masimanimba", "Low_Prevalence_Zone", "High_Prevalence_Zone"))

#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

my_comparisons <- list(c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone")) #KO5350
#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001

t <- ggplot(K.tr.DF,aes(x = Status,y = K05350)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab("rel. abund. of K05350: beta-glucosidase [EC:3.2.1.21]") + stat_boxplot(geom ='errorbar')
t <- t + geom_jitter(position=position_jitter(0.2), size = 0.2)
t <- t + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = konzo_color) + theme(plot.title = element_blank(), legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 6), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 4.5), axis.title.x = element_blank())
t <- t + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)


tiff(filename = "Kinshasa_Konzo3_KO_K05350_betaglucosidase_BoxPlot.tiff", width = 3.5, height = 3.5, units = "in", res = 600)
t
dev.off()

###
tiff(filename = "Kinshasa_Konzo3_LAB_K05350_BoxPlot.tiff", width = 6, height = 4, units = "in", res = 600)
ggarrange(lab,t,labels = c("A","B"), widths = c(3.5, 2.5), ncol = 2, nrow = 1, font.label = list(size = 7))
dev.off()
                                    
                                    
#E. coli
ec <- ggplot(S.tr.DF,aes(x = Status, y = S.tr.DF$Escherichia.coli)) + 
    geom_boxplot(aes(fill = Status), outlier.size = 0.2, fatten = 0.5) + theme_classic() + ylab(expression(paste("rel. abund. of ", italic("Escherichia coli"))))
ec <- ec + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = konzo_color) + theme(plot.title = element_blank(), legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 6), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 7), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 6), axis.title.x = element_blank())
ec <- ec + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)
                                    

#Rhodanese (K01011)
#thiosulfate/3-mercaptopyruvate sulfurtransferase [EC:2.8.1.1, 2.8.1.2]                                    
                                    
K <- KonzoData_KO_tr

K.tr_META <- as.data.frame(K@sam_data)
K.tr_OTU <- as.data.frame(t(K@otu_table))
K.tr.DF <- cbind(K.tr_OTU, K.tr_META$Status)
K.tr.DF <- cbind(K.tr.DF, K.tr_META$Geography)

colnames(K.tr.DF)[colnames(K.tr.DF)=="K.tr_META$Status"] <- "Status"
colnames(K.tr.DF)[colnames(K.tr.DF)=="K.tr_META$Geography"] <- "Geography"

for (i in nrow(K.tr.DF))
{K.tr.DF[i,]$Status <- KonzoData_KO_tr@sam_data[rownames(K.tr.DF[i,]),]$Status
}
K.tr.DF$Status <- factor(K.tr.DF$Status, levels = c("Kinshasa", "Masimanimba", "Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

for (i in nrow(K.tr.DF))
{K.tr.DF[i,]$Geography <- KonzoData_KO_tr@sam_data[rownames(K.tr.DF[i,]),]$Geography
}
K.tr.DF$Geography <- factor(K.tr.DF$Geography, levels = c("Kinshasa", "Masimanimba", "Low_Prevalence_Zone", "High_Prevalence_Zone"))
                                    
                                    
pw_wt <-  pairwise.wilcox.test(K.tr.DF$K01011, K.tr.DF$Status, p.adj = "none")
#capture.output(pw_wt, file="KinshasaKonzo3_pairwiseMWW_K01011_rhodanese.txt")                                                
                                    
                                    
#my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Kinshasa", "Unaffected_Low_Prevalence_Zone"), c("Kinshasa", "Konzo_Low_Prevalence_Zone"), c("Kinshasa", "Unaffected_High_Prevalence_Zone"), c("Kinshasa", "Konzo_High_Prevalence_Zone"), 
                        #c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone"), 
                       #c("Unaffected_Low_Prevalence_Zone", "Unaffected_High_Prevalence_Zone"), c("Konzo_Low_Prevalence_Zone", "Konzo_High_Prevalence_Zone"), 
                        #c("Unaffected_Low_Prevalence_Zone", "Konzo_Low_Prevalence_Zone"), c("Unaffected_High_Prevalence_Zone", "Konzo_High_Prevalence_Zone"))

my_comparisons <- list( c("Kinshasa", "Masimanimba"), c("Masimanimba", "Unaffected_Low_Prevalence_Zone"), c("Masimanimba", "Konzo_Low_Prevalence_Zone"), c("Masimanimba", "Unaffected_High_Prevalence_Zone"), c("Masimanimba", "Konzo_High_Prevalence_Zone")) #K01011

#ns: p > 0.05
#*: p <= 0.05
#**: p <= 0.01
#***: p <= 0.001
#****: p <= 0.0001

r <- ggplot(K.tr.DF,aes(x = Status,y = K01011)) + 
    geom_boxplot(aes(fill = Status),outlier.shape = NA, fatten = 0.5) + theme_classic() + ylab(expression(paste("rel. abund. of K01011: \nthiosulfate/3-mercaptopyruvate sulfurtransferase \n[EC:2.8.1.1, 2.8.1.2]"))) + stat_boxplot(geom ='errorbar')
r <- r + geom_jitter(position=position_jitter(0.2), size = 0.3)
r <- r + theme(legend.position="NA") + scale_x_discrete(labels= SSSL) + scale_fill_manual(values = konzo_color) + theme(plot.title = element_blank(), legend.key.size = unit(.4, "cm"), legend.text = element_text(size = 6), legend.title = element_blank()) + 
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6), axis.title.y = element_text(size = 4), axis.title.x = element_blank())
r <- r + stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test", size = 2)

                                    
tiff(filename = "Kinshasa_Konzo3_Ecoli_K01011_BoxPlot.tiff", width = 6, height = 4, units = "in", res = 600)
ggarrange(ec,r,labels = c("C","D"), ncol = 2, nrow = 1, font.label = list(size = 7))
dev.off() 
                                    
                                    
#tiff(filename = "Kinshasa_Konzo3_Lab_Ecoli_Functional_BoxPlot.tiff", width = 7, height = 5, units = "in", res = 600)
#ggarrange(lab, ec, t ,r,labels = c("A", "C","B","D"), heights = c(3,2), ncol = 2, nrow = 2, font.label = list(size = 7))
#dev.off()    
                                    
tiff(filename = "Kinshasa_Konzo3_Lab_Functional_BoxPlot.tiff", width = 7, height = 3.5, units = "in", res = 600)
ggarrange(lab, t , r,labels = c("A","B","C"), widths = c(1.35,1,1), ncol = 3, nrow = 1, font.label = list(size = 7))
dev.off()  
                            
                            
                            
                            
                            
###THE END ######
                            
                            
                                                                         
                                                  
