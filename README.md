# Gut_Microbiome_Konzo
The Gut Microbiome in Konzo

The code outlined here assesses the bacterial profile of the gut microbiome of individuals from the Democratic Republic of Congo.

Note: The shell scripts are there for reference and will need to be adjusted depending on input folder. However the parameters are the same for all processed data. 

Note: BMTagger, skewer, Kraken2, and Bracken were run on linux based OS, while R was run on a macOS with 32GB memory.

Note: The code is commented up front for clarification, but most of the code is repetitive with changes depending on input. Comments are added whenever explaination may be necessary and refer back when needed. Additionally, phylum input files are available for reference to have a better idea of what the input looks like as the files are separate and have more info than necessary/provided in Supplemental files in the paper submission. The Supplemental file 1 includes all information needed to generate the same data output. 

1. Run bmtagger.sh (version 1.1.0) to remove human reads by aligning to the human genome reference, hg38. The bmtagger.conf file may be necessary even if you have installed all dependencies. To successfully run BMTagger required programs, about 8.5Gb memory and around three times harddisk space to store index data is needed. Follow installation and running guidelines here: https://hmpdacc.org/hmp/doc/HumanSequenceRemoval_SOP.pdf  
2. Run skewer.sh (version 0.2.2) on the fastq files to remove Illumina sequencing adapters. 
3. Run kraken2.sh on post skewer using Kraken2 (version 2.0.8-beta) standard genome libraries (kraken2-build.sh; built on March 3rd, 2019), which includes human, viral, bacteria, and archae. Construction of the kraken database can require about 100GB of space, and enough free memory is needed in RAM to build the databass. Thus for ~29GB of standard database, you need a little more RAM than that to accomodate. The detailed manual for the program installation and use can be found here: https://github.com/DerrickWood/kraken2/wiki/Manual
4. Run bracken.sh on the classified reads from Kraken2 for secondary analysis. The manual for installation and use can be found here: https://ccb.jhu.edu/software/bracken/index.shtml?t=manual
5. Extract the Bracken files that are in Kraken report format (should be in the input folder) for downstream analysis in R since they are compatible with pavian. 
6. Use pavain to consolidate the report files for all samples into one large file (will need to converted to csv from tsv). Additionally, using the pavian filter feature, get different taxon (phylum, class, order, family, genus, species) for each of the different domains (here only Bacteria are analyzed). The lineage column is removed and the names of samples are changed to they match those in the Konzo_meta file prior to analysis in R. Additionaly information on pavian can be found here: https://ccb.jhu.edu/software/pavian/index.shtml
7. Install any additional necessary packages as listed in Code.R (R version 3.6.3), and follow the script for further analysis.
8. Related figure making is also included in the code for reference. Gimp was also used to combine some figures into a final layout. 


