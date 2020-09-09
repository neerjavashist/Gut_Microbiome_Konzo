# Gut_Microbiome_Konzo
The Gut Microbiome in Konzo

The code outlined here assesses the bacterial profile of the gut microbiome of individuals from the Democratic Republic of Congo.
The shell scripts are there for referece and will need to be adjusted depending on input folder. However the parameters are the same for all processed data 

1. Run bmtagger.sh to remove human reads using hg38 and the needed bmtagger.conf file. 
2. Run skewer.sh on the fastq files to remove Illumina sequencing adapters. 
3. Run kraken2.sh on the fastq post skewer using Kraken2 standard genome libraries (kraken2-build.sh; built on March 3rd, 2019), which includes human, viral, bacteria, and archae. The detailed manual for the program installation and use can be found here: https://github.com/DerrickWood/kraken2/wiki/Manual
4. Run bracken.sh on the classified reads from Kraken2 for secondary analysis and the manual for installation and use can be found here: https://ccb.jhu.edu/software/bracken/index.shtml?t=manual
5. Extract the Bracken files that in Kraken report format (will be in the input folder) for downstream analysis in R since they are compatible with pavian. 
6. Use pavain to consolidate the report files for all samples into one large file. Additionally, using the pavian filter feature, get different taxon (phylum, class, order, family, genus, species) for each of the different domains (here only Bacteria are analyzed). The lineage column is removed and the names of samples are changed to they match those in the Konzo_meta file prior to analysis in R. Additionaly information on pavian can be found here: https://ccb.jhu.edu/software/pavian/index.shtml
7. Install any additional necessary packages as listed in Code.R, and follow the script for further analysis.
8. Figure making is also included in the code. However, Gimp was used to combine some figures into a final layout. 

Note: The code is commented up front for clarification, but most of the code is repetitive with changes depedning on input. Additionally, phylum input files are available for reference to have a better idea of what the input looks like as the files are separate and have more info than necessary/provided in Supplemental files in the paper submission. 
