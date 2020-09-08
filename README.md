# Gut_Microbiome_Konzo
The Gut Microbiome in Konzo

The code outlined here assesses the bacterial profile of the gut microbiome of individuals from the Democratic Republic of Congo.

1. Run bmtagger.sh to remove human reads using hg38 and the needed bmtagger.conf file. 
2. Run skewer.sh on the fastq files to remove Illumina sequencing adapters. 
3. Run kraken2.sh on the fastq post skewer using Kraken2 standard genome libraries, which include human, viral, bacteria, and archae. The detailed manual for the program can be found here: https://github.com/DerrickWood/kraken2/wiki/Manual
4. Run bracken.sh on the classified reads from Kraken2 for secondary analysis.
5. Extract the Braken files that in Kraken report format (will be in the input folder) for downstream analysis in R since they are compatible with pavian. 
6. Use pavain to consolidate the report files for all samples into one large file. Additionally, using the pavian filter feature, get different taxon (phylum, class, order, family, genus, species) for each of the different domains (here only Bacteria are analyzed). 
7. Install any additional necessary packages as listed in Code.R, and follow the script for further analysis. 
