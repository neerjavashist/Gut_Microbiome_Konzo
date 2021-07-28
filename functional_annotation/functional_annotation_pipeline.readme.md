# Functional annotation for konzo metagenomic analysis
## Tools required
- FASTP
***https://github.com/OpenGene/fastp***

- KRAKEN2
***https://github.com/DerrickWood/kraken2***

- BOWTIE2

***http://bowtie-bio.sourceforge.net/bowtie2/index.shtml***

```
cd Gut_Microbiome_Konzo/functional_annotation
```

# Download genomes for database creation

***https://github.com/kblin/ncbi-genome-download***

```
cd db
ncbi-genome-download --format fasta,assembly-report --assembly-level chromosome  --parallel 4  bacteria
ncbi-genome-download --format fasta,assembly-report --assembly-level chromosome  --parallel 2 viral
ncbi-genome-download --format fasta,assembly-report --assembly-level chromosome  --parallel 2 fungi
ncbi-genome-download --format fasta,assembly-report --assembly-level chromosome  --parallel 2  archaea
ncbi-genome-download --format fasta,assembly-report --taxid 9606 --assembly-level chromosome vertebrate_mammalian
ncbi-genome-download --format fasta,assembly-report --assembly-level complete  archaea
ncbi-genome-download --format fasta,assembly-report --assembly-level complete  bacteria
ncbi-genome-download  --format fasta,assembly-report --assembly-level complete  --parallel 4  viral
ncbi-genome-download  --format fasta,assembly-report --assembly-level complete  --parallel 2  fungi
ncbi-genome-download --format fasta,assembly-report --assembly-level scaffold  --parallel 2  bacteria
ncbi-genome-download --format fasta,assembly-report --assembly-level scaffold  --parallel 2 archaea
ncbi-genome-download --format fasta,assembly-report --assembly-level scaffold  --parallel 2 viral
ncbi-genome-download --format fasta,assembly-report --assembly-level scaffold  --parallel 2 fungi
```

# Taxonomic annotation with Database creation kraken2 kraken2/2.0.8-beta
```
for f in */refseq/*/*/*.fna.gz; do zcat "$f" >> library/gall_genomes.DB.fna
```
the entire fasta file has a size of approximately 330 Go to process the data you need to have 330 minimum RAM

```
cd ..
kraken2-build --standard --db db/
```


# bowtie2 db creation of kegg microbial gene database
***https://www.genome.jp/kegg/***
You have to buy the kegg db in order to obtain the following data:

- ko_genes.list
- ko_module.list
- ko_pathway.list
- kegg microbial gene sequences file

once you have the required fasta genes of all microbial genes:
```
cd db/
bowtie2-build --large-index --threads 40 {fasta_genes} bowtie2_keggDB
cd ..
```

# Quality trimming with fastp/0.20.0

```
fastp -i rawdata/L0050-34-22-38-Vilain-180822_S34_R1_001.fastq -I rawdata/L0050-34-22-38-Vilain-180822_S34_R2_001.fastq -o preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R1_001.good.fq -O preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R2_001.good.fq -w 2
fastp -i rawdata/L0070-54-10-079-Vilain-180822_S54_R1_001.fastq -I rawdata/L0070-54-10-079-Vilain-180822_S54_R2_001.fastq -o preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R1_001.good.fq -O preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R2_001.good.fq -w 2
fastp -i rawdata/L0081-65-10-079-Vilain-180822_S65_R1_001.fastq -I rawdata/L0081-65-10-079-Vilain-180822_S65_R2_001.fastq -o preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R1_001.good.fq -O preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R2_001.good.fq -w 2
fastp -i rawdata/Konzo10_S91_ME_L001_R1_001.fastq -I rawdata/Konzo10_S91_ME_L001_R2_001.fastq -o preprocessing/good/Konzo10_S91_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo10_S91_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo11_S90_ME_L001_R1_001.fastq -I rawdata/Konzo11_S90_ME_L001_R2_001.fastq -o preprocessing/good/Konzo11_S90_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo11_S90_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo12_S89_ME_L001_R1_001.fastq -I rawdata/Konzo12_S89_ME_L001_R2_001.fastq -o preprocessing/good/Konzo12_S89_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo12_S89_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo13_S88_ME_L001_R1_001.fastq -I rawdata/Konzo13_S88_ME_L001_R2_001.fastq -o preprocessing/good/Konzo13_S88_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo13_S88_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo14_S2_ME_L001_R1_001.fastq -I rawdata/Konzo14_S2_ME_L001_R2_001.fastq -o preprocessing/good/Konzo14_S2_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo14_S2_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo15_S3_ME_L001_R1_001.fastq -I rawdata/Konzo15_S3_ME_L001_R2_001.fastq -o preprocessing/good/Konzo15_S3_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo15_S3_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo16_S4_ME_L001_R1_001.fastq -I rawdata/Konzo16_S4_ME_L001_R2_001.fastq -o preprocessing/good/Konzo16_S4_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo16_S4_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo17_S5_ME_L001_R1_001.fastq -I rawdata/Konzo17_S5_ME_L001_R2_001.fastq -o preprocessing/good/Konzo17_S5_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo17_S5_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo18_S87_ME_L001_R1_001.fastq -I rawdata/Konzo18_S87_ME_L001_R2_001.fastq -o preprocessing/good/Konzo18_S87_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo18_S87_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo19_S86_ME_L001_R1_001.fastq -I rawdata/Konzo19_S86_ME_L001_R2_001.fastq -o preprocessing/good/Konzo19_S86_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo19_S86_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo1_S1_ME_L001_R1_001.fastq -I rawdata/Konzo1_S1_ME_L001_R2_001.fastq -o preprocessing/good/Konzo1_S1_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo1_S1_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo20_S85_ME_L001_R1_001.fastq -I rawdata/Konzo20_S85_ME_L001_R2_001.fastq -o preprocessing/good/Konzo20_S85_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo20_S85_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo21_S84_ME_L001_R1_001.fastq -I rawdata/Konzo21_S84_ME_L001_R2_001.fastq -o preprocessing/good/Konzo21_S84_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo21_S84_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo22_S83_ME_L001_R1_001.fastq -I rawdata/Konzo22_S83_ME_L001_R2_001.fastq -o preprocessing/good/Konzo22_S83_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo22_S83_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo23_S82_ME_L001_R1_001.fastq -I rawdata/Konzo23_S82_ME_L001_R2_001.fastq -o preprocessing/good/Konzo23_S82_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo23_S82_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo24_S81_ME_L001_R1_001.fastq -I rawdata/Konzo24_S81_ME_L001_R2_001.fastq -o preprocessing/good/Konzo24_S81_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo24_S81_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo25_S80_ME_L001_R1_001.fastq -I rawdata/Konzo25_S80_ME_L001_R2_001.fastq -o preprocessing/good/Konzo25_S80_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo25_S80_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo26_S7_ME_L001_R1_001.fastq -I rawdata/Konzo26_S7_ME_L001_R2_001.fastq -o preprocessing/good/Konzo26_S7_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo26_S7_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo27_S8_ME_L001_R1_001.fastq -I rawdata/Konzo27_S8_ME_L001_R2_001.fastq -o preprocessing/good/Konzo27_S8_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo27_S8_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo28_S9_ME_L001_R1_001.fastq -I rawdata/Konzo28_S9_ME_L001_R2_001.fastq -o preprocessing/good/Konzo28_S9_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo28_S9_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo29_S10_ME_L001_R1_001.fastq -I rawdata/Konzo29_S10_ME_L001_R2_001.fastq -o preprocessing/good/Konzo29_S10_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo29_S10_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo2_S6_ME_L001_R1_001.fastq -I rawdata/Konzo2_S6_ME_L001_R2_001.fastq -o preprocessing/good/Konzo2_S6_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo2_S6_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo30_S79_ME_L001_R1_001.fastq -I rawdata/Konzo30_S79_ME_L001_R2_001.fastq -o preprocessing/good/Konzo30_S79_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo30_S79_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo31_S78_ME_L001_R1_001.fastq -I rawdata/Konzo31_S78_ME_L001_R2_001.fastq -o preprocessing/good/Konzo31_S78_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo31_S78_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo32_S12_ME_L001_R1_001.fastq -I rawdata/Konzo32_S12_ME_L001_R2_001.fastq -o preprocessing/good/Konzo32_S12_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo32_S12_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo33_S77_ME_L001_R1_001.fastq -I rawdata/Konzo33_S77_ME_L001_R2_001.fastq -o preprocessing/good/Konzo33_S77_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo33_S77_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo34_S76_ME_L001_R1_001.fastq -I rawdata/Konzo34_S76_ME_L001_R2_001.fastq -o preprocessing/good/Konzo34_S76_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo34_S76_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo35_S75_ME_L001_R1_001.fastq -I rawdata/Konzo35_S75_ME_L001_R2_001.fastq -o preprocessing/good/Konzo35_S75_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo35_S75_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo36_S74_ME_L001_R1_001.fastq -I rawdata/Konzo36_S74_ME_L001_R2_001.fastq -o preprocessing/good/Konzo36_S74_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo36_S74_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo37_S73_ME_L001_R1_001.fastq -I rawdata/Konzo37_S73_ME_L001_R2_001.fastq -o preprocessing/good/Konzo37_S73_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo37_S73_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/Konzo38_S72_ME_L001_R1_001.fastq -I rawdata/Konzo38_S72_ME_L001_R2_001.fastq -o preprocessing/good/Konzo38_S72_ME_L001_R1_001.good.fq -O preprocessing/good/Konzo38_S72_ME_L001_R2_001.good.fq -w 2
fastp -i rawdata/L0049-33-10-76-Vilain-180822_S33_R1_001.fastq -I rawdata/L0049-33-10-76-Vilain-180822_S33_R2_001.fastq -o preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R1_001.good.fq -O preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R2_001.good.fq -w 2
fastp -i rawdata/L0060-44-22-38-Vilain-180822_S44_R1_001.fastq -I rawdata/L0060-44-22-38-Vilain-180822_S44_R2_001.fastq -o preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R1_001.good.fq -O preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R2_001.good.fq -w 2
fastp -i rawdata/L0063-47-10-81-Vilain-180822_S47_R1_001.fastq -I rawdata/L0063-47-10-81-Vilain-180822_S47_R2_001.fastq -o preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R1_001.good.fq -O preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R2_001.good.fq -w 2
fastp -i rawdata/L0067-51-10-76-Vilain-180822_S51_R1_001.fastq -I rawdata/L0067-51-10-76-Vilain-180822_S51_R2_001.fastq -o preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R1_001.good.fq -O preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R2_001.good.fq -w 2
fastp -i rawdata/L337-A34-Vilain-181127_S32_R1_001.fastq -I rawdata/L337-A34-Vilain-181127_S32_R2_001.fastq -o preprocessing/good/L337-A34-Vilain-181127_S32_R1_001.good.fq -O preprocessing/good/L337-A34-Vilain-181127_S32_R2_001.good.fq -w 2
fastp -i rawdata/L338-A20-Vilain-181127_S33_R1_001.fastq -I rawdata/L338-A20-Vilain-181127_S33_R2_001.fastq -o preprocessing/good/L338-A20-Vilain-181127_S33_R1_001.good.fq -O preprocessing/good/L338-A20-Vilain-181127_S33_R2_001.good.fq -w 2
fastp -i rawdata/L339-A31-Vilain-181127_S34_R1_001.fastq -I rawdata/L339-A31-Vilain-181127_S34_R2_001.fastq -o preprocessing/good/L339-A31-Vilain-181127_S34_R1_001.good.fq -O preprocessing/good/L339-A31-Vilain-181127_S34_R2_001.good.fq -w 2
fastp -i rawdata/L340-A33-Vilain-181127_S35_R1_001.fastq -I rawdata/L340-A33-Vilain-181127_S35_R2_001.fastq -o preprocessing/good/L340-A33-Vilain-181127_S35_R1_001.good.fq -O preprocessing/good/L340-A33-Vilain-181127_S35_R2_001.good.fq -w 2
fastp -i rawdata/L341-A26-Vilain-181127_S36_R1_001.fastq -I rawdata/L341-A26-Vilain-181127_S36_R2_001.fastq -o preprocessing/good/L341-A26-Vilain-181127_S36_R1_001.good.fq -O preprocessing/good/L341-A26-Vilain-181127_S36_R2_001.good.fq -w 2
fastp -i rawdata/L342-A13-Vilain-181127_S37_R1_001.fastq -I rawdata/L342-A13-Vilain-181127_S37_R2_001.fastq -o preprocessing/good/L342-A13-Vilain-181127_S37_R1_001.good.fq -O preprocessing/good/L342-A13-Vilain-181127_S37_R2_001.good.fq -w 2
fastp -i rawdata/L343-A32-Vilain-181127_S38_R1_001.fastq -I rawdata/L343-A32-Vilain-181127_S38_R2_001.fastq -o preprocessing/good/L343-A32-Vilain-181127_S38_R1_001.good.fq -O preprocessing/good/L343-A32-Vilain-181127_S38_R2_001.good.fq -w 2
fastp -i rawdata/L344-A3-Vilain-181127_S39_R1_001.fastq -I rawdata/L344-A3-Vilain-181127_S39_R2_001.fastq -o preprocessing/good/L344-A3-Vilain-181127_S39_R1_001.good.fq -O preprocessing/good/L344-A3-Vilain-181127_S39_R2_001.good.fq -w 2
fastp -i rawdata/L345-A6-Vilain-181127_S40_R1_001.fastq -I rawdata/L345-A6-Vilain-181127_S40_R2_001.fastq -o preprocessing/good/L345-A6-Vilain-181127_S40_R1_001.good.fq -O preprocessing/good/L345-A6-Vilain-181127_S40_R2_001.good.fq -w 2
```

# taxonomic annotation with  kraken2/2.0.8-beta

```
kraken2 -t 40 -db db/ --output output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.bigDBkraken2.txt --report output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.bigDBkraken2_report.tsv --paired preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R1_001.good.fq preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.bigDBkraken2.txt --report output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.bigDBkraken2_report.tsv --paired preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R1_001.good.fq preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.bigDBkraken2.txt --report output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.bigDBkraken2_report.tsv --paired preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R1_001.good.fq preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo10_S91_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo10_S91_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo10_S91_ME_L001_R1_001.good.fq preprocessing/good/Konzo10_S91_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo11_S90_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo11_S90_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo11_S90_ME_L001_R1_001.good.fq preprocessing/good/Konzo11_S90_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo12_S89_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo12_S89_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo12_S89_ME_L001_R1_001.good.fq preprocessing/good/Konzo12_S89_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo13_S88_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo13_S88_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo13_S88_ME_L001_R1_001.good.fq preprocessing/good/Konzo13_S88_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo14_S2_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo14_S2_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo14_S2_ME_L001_R1_001.good.fq preprocessing/good/Konzo14_S2_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo15_S3_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo15_S3_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo15_S3_ME_L001_R1_001.good.fq preprocessing/good/Konzo15_S3_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo16_S4_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo16_S4_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo16_S4_ME_L001_R1_001.good.fq preprocessing/good/Konzo16_S4_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo17_S5_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo17_S5_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo17_S5_ME_L001_R1_001.good.fq preprocessing/good/Konzo17_S5_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo18_S87_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo18_S87_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo18_S87_ME_L001_R1_001.good.fq preprocessing/good/Konzo18_S87_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo19_S86_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo19_S86_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo19_S86_ME_L001_R1_001.good.fq preprocessing/good/Konzo19_S86_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo1_S1_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo1_S1_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo1_S1_ME_L001_R1_001.good.fq preprocessing/good/Konzo1_S1_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo20_S85_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo20_S85_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo20_S85_ME_L001_R1_001.good.fq preprocessing/good/Konzo20_S85_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo21_S84_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo21_S84_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo21_S84_ME_L001_R1_001.good.fq preprocessing/good/Konzo21_S84_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo22_S83_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo22_S83_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo22_S83_ME_L001_R1_001.good.fq preprocessing/good/Konzo22_S83_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo23_S82_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo23_S82_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo23_S82_ME_L001_R1_001.good.fq preprocessing/good/Konzo23_S82_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo24_S81_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo24_S81_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo24_S81_ME_L001_R1_001.good.fq preprocessing/good/Konzo24_S81_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo25_S80_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo25_S80_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo25_S80_ME_L001_R1_001.good.fq preprocessing/good/Konzo25_S80_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo26_S7_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo26_S7_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo26_S7_ME_L001_R1_001.good.fq preprocessing/good/Konzo26_S7_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo27_S8_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo27_S8_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo27_S8_ME_L001_R1_001.good.fq preprocessing/good/Konzo27_S8_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo28_S9_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo28_S9_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo28_S9_ME_L001_R1_001.good.fq preprocessing/good/Konzo28_S9_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo29_S10_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo29_S10_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo29_S10_ME_L001_R1_001.good.fq preprocessing/good/Konzo29_S10_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo2_S6_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo2_S6_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo2_S6_ME_L001_R1_001.good.fq preprocessing/good/Konzo2_S6_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo30_S79_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo30_S79_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo30_S79_ME_L001_R1_001.good.fq preprocessing/good/Konzo30_S79_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo31_S78_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo31_S78_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo31_S78_ME_L001_R1_001.good.fq preprocessing/good/Konzo31_S78_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo32_S12_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo32_S12_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo32_S12_ME_L001_R1_001.good.fq preprocessing/good/Konzo32_S12_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo33_S77_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo33_S77_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo33_S77_ME_L001_R1_001.good.fq preprocessing/good/Konzo33_S77_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo34_S76_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo34_S76_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo34_S76_ME_L001_R1_001.good.fq preprocessing/good/Konzo34_S76_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo35_S75_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo35_S75_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo35_S75_ME_L001_R1_001.good.fq preprocessing/good/Konzo35_S75_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo36_S74_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo36_S74_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo36_S74_ME_L001_R1_001.good.fq preprocessing/good/Konzo36_S74_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo37_S73_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo37_S73_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo37_S73_ME_L001_R1_001.good.fq preprocessing/good/Konzo37_S73_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/Konzo38_S72_ME_L001.bigDBkraken2.txt --report output/kraken2_output/Konzo38_S72_ME_L001.bigDBkraken2_report.tsv --paired preprocessing/good/Konzo38_S72_ME_L001_R1_001.good.fq preprocessing/good/Konzo38_S72_ME_L001_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.bigDBkraken2.txt --report output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.bigDBkraken2_report.tsv --paired preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R1_001.good.fq preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.bigDBkraken2.txt --report output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.bigDBkraken2_report.tsv --paired preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R1_001.good.fq preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.bigDBkraken2.txt --report output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.bigDBkraken2_report.tsv --paired preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R1_001.good.fq preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.bigDBkraken2.txt --report output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.bigDBkraken2_report.tsv --paired preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R1_001.good.fq preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L337-A34-Vilain-181127_S32.bigDBkraken2.txt --report output/kraken2_output/L337-A34-Vilain-181127_S32.bigDBkraken2_report.tsv --paired preprocessing/good/L337-A34-Vilain-181127_S32_R1_001.good.fq preprocessing/good/L337-A34-Vilain-181127_S32_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L338-A20-Vilain-181127_S33.bigDBkraken2.txt --report output/kraken2_output/L338-A20-Vilain-181127_S33.bigDBkraken2_report.tsv --paired preprocessing/good/L338-A20-Vilain-181127_S33_R1_001.good.fq preprocessing/good/L338-A20-Vilain-181127_S33_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L339-A31-Vilain-181127_S34.bigDBkraken2.txt --report output/kraken2_output/L339-A31-Vilain-181127_S34.bigDBkraken2_report.tsv --paired preprocessing/good/L339-A31-Vilain-181127_S34_R1_001.good.fq preprocessing/good/L339-A31-Vilain-181127_S34_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L340-A33-Vilain-181127_S35.bigDBkraken2.txt --report output/kraken2_output/L340-A33-Vilain-181127_S35.bigDBkraken2_report.tsv --paired preprocessing/good/L340-A33-Vilain-181127_S35_R1_001.good.fq preprocessing/good/L340-A33-Vilain-181127_S35_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L341-A26-Vilain-181127_S36.bigDBkraken2.txt --report output/kraken2_output/L341-A26-Vilain-181127_S36.bigDBkraken2_report.tsv --paired preprocessing/good/L341-A26-Vilain-181127_S36_R1_001.good.fq preprocessing/good/L341-A26-Vilain-181127_S36_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L342-A13-Vilain-181127_S37.bigDBkraken2.txt --report output/kraken2_output/L342-A13-Vilain-181127_S37.bigDBkraken2_report.tsv --paired preprocessing/good/L342-A13-Vilain-181127_S37_R1_001.good.fq preprocessing/good/L342-A13-Vilain-181127_S37_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L343-A32-Vilain-181127_S38.bigDBkraken2.txt --report output/kraken2_output/L343-A32-Vilain-181127_S38.bigDBkraken2_report.tsv --paired preprocessing/good/L343-A32-Vilain-181127_S38_R1_001.good.fq preprocessing/good/L343-A32-Vilain-181127_S38_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L344-A3-Vilain-181127_S39.bigDBkraken2.txt --report output/kraken2_output/L344-A3-Vilain-181127_S39.bigDBkraken2_report.tsv --paired preprocessing/good/L344-A3-Vilain-181127_S39_R1_001.good.fq preprocessing/good/L344-A3-Vilain-181127_S39_R2_001.good.fq
kraken2 -t 40 -db db/ --output output/kraken2_output/L345-A6-Vilain-181127_S40.bigDBkraken2.txt --report output/kraken2_output/L345-A6-Vilain-181127_S40.bigDBkraken2_report.tsv --paired preprocessing/good/L345-A6-Vilain-181127_S40_R1_001.good.fq preprocessing/good/L345-A6-Vilain-181127_S40_R2_001.good.fq
```

# Sort read ID affiliated to Homo sapiens

```
 grep '^C' output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.non_human_reads.txt
 grep '^C' output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.non_human_reads.txt
 grep '^C' output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo10_S91_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo10_S91_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo11_S90_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo11_S90_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo12_S89_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo12_S89_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo13_S88_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo13_S88_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo14_S2_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo14_S2_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo15_S3_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo15_S3_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo16_S4_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo16_S4_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo17_S5_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo17_S5_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo18_S87_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo18_S87_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo19_S86_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo19_S86_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo1_S1_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo1_S1_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo20_S85_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo20_S85_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo21_S84_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo21_S84_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo22_S83_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo22_S83_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo23_S82_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo23_S82_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo24_S81_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo24_S81_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo25_S80_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo25_S80_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo26_S7_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo26_S7_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo27_S8_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo27_S8_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo28_S9_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo28_S9_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo29_S10_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo29_S10_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo2_S6_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo2_S6_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo30_S79_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo30_S79_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo31_S78_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo31_S78_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo32_S12_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo32_S12_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo33_S77_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo33_S77_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo34_S76_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo34_S76_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo35_S75_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo35_S75_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo36_S74_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo36_S74_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo37_S73_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo37_S73_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/Konzo38_S72_ME_L001.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/Konzo38_S72_ME_L001.non_human_reads.txt
 grep '^C' output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.non_human_reads.txt
 grep '^C' output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.non_human_reads.txt
 grep '^C' output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.non_human_reads.txt
 grep '^C' output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.non_human_reads.txt
 grep '^C' output/kraken2_output/L337-A34-Vilain-181127_S32.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L337-A34-Vilain-181127_S32.non_human_reads.txt
 grep '^C' output/kraken2_output/L338-A20-Vilain-181127_S33.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L338-A20-Vilain-181127_S33.non_human_reads.txt
 grep '^C' output/kraken2_output/L339-A31-Vilain-181127_S34.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L339-A31-Vilain-181127_S34.non_human_reads.txt
 grep '^C' output/kraken2_output/L340-A33-Vilain-181127_S35.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L340-A33-Vilain-181127_S35.non_human_reads.txt
 grep '^C' output/kraken2_output/L341-A26-Vilain-181127_S36.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L341-A26-Vilain-181127_S36.non_human_reads.txt
 grep '^C' output/kraken2_output/L342-A13-Vilain-181127_S37.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L342-A13-Vilain-181127_S37.non_human_reads.txt
 grep '^C' output/kraken2_output/L343-A32-Vilain-181127_S38.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L343-A32-Vilain-181127_S38.non_human_reads.txt
 grep '^C' output/kraken2_output/L344-A3-Vilain-181127_S39.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L344-A3-Vilain-181127_S39.non_human_reads.txt
 grep '^C' output/kraken2_output/L345-A6-Vilain-181127_S40.bigDBkraken2.txt |grep -v '	9606 |cut -f 2 > output/kraken2_output/L345-A6-Vilain-181127_S40.non_human_reads.txt
```

# filter fastq to recrute non human reads

```
perl scripts/filter_fastq.pl -R=preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R1_001.good.fq -out=output/function/L0050-34-22-38-Vilain-180822_S34_R1.human_trimmed.fq -human=output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R1_001.good.fq -out=output/function/L0070-54-10-079-Vilain-180822_S54_R1.human_trimmed.fq -human=output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R1_001.good.fq -out=output/function/L0081-65-10-079-Vilain-180822_S65_R1.human_trimmed.fq -human=output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo10_S91_ME_L001_R1_001.good.fq -out=output/function/Konzo10_S91_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo10_S91_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo11_S90_ME_L001_R1_001.good.fq -out=output/function/Konzo11_S90_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo11_S90_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo12_S89_ME_L001_R1_001.good.fq -out=output/function/Konzo12_S89_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo12_S89_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo13_S88_ME_L001_R1_001.good.fq -out=output/function/Konzo13_S88_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo13_S88_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo14_S2_ME_L001_R1_001.good.fq -out=output/function/Konzo14_S2_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo14_S2_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo15_S3_ME_L001_R1_001.good.fq -out=output/function/Konzo15_S3_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo15_S3_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo16_S4_ME_L001_R1_001.good.fq -out=output/function/Konzo16_S4_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo16_S4_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo17_S5_ME_L001_R1_001.good.fq -out=output/function/Konzo17_S5_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo17_S5_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo18_S87_ME_L001_R1_001.good.fq -out=output/function/Konzo18_S87_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo18_S87_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo19_S86_ME_L001_R1_001.good.fq -out=output/function/Konzo19_S86_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo19_S86_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo1_S1_ME_L001_R1_001.good.fq -out=output/function/Konzo1_S1_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo1_S1_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo20_S85_ME_L001_R1_001.good.fq -out=output/function/Konzo20_S85_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo20_S85_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo21_S84_ME_L001_R1_001.good.fq -out=output/function/Konzo21_S84_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo21_S84_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo22_S83_ME_L001_R1_001.good.fq -out=output/function/Konzo22_S83_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo22_S83_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo23_S82_ME_L001_R1_001.good.fq -out=output/function/Konzo23_S82_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo23_S82_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo24_S81_ME_L001_R1_001.good.fq -out=output/function/Konzo24_S81_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo24_S81_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo25_S80_ME_L001_R1_001.good.fq -out=output/function/Konzo25_S80_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo25_S80_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo26_S7_ME_L001_R1_001.good.fq -out=output/function/Konzo26_S7_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo26_S7_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo27_S8_ME_L001_R1_001.good.fq -out=output/function/Konzo27_S8_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo27_S8_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo28_S9_ME_L001_R1_001.good.fq -out=output/function/Konzo28_S9_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo28_S9_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo29_S10_ME_L001_R1_001.good.fq -out=output/function/Konzo29_S10_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo29_S10_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo2_S6_ME_L001_R1_001.good.fq -out=output/function/Konzo2_S6_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo2_S6_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo30_S79_ME_L001_R1_001.good.fq -out=output/function/Konzo30_S79_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo30_S79_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo31_S78_ME_L001_R1_001.good.fq -out=output/function/Konzo31_S78_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo31_S78_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo32_S12_ME_L001_R1_001.good.fq -out=output/function/Konzo32_S12_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo32_S12_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo33_S77_ME_L001_R1_001.good.fq -out=output/function/Konzo33_S77_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo33_S77_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo34_S76_ME_L001_R1_001.good.fq -out=output/function/Konzo34_S76_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo34_S76_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo35_S75_ME_L001_R1_001.good.fq -out=output/function/Konzo35_S75_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo35_S75_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo36_S74_ME_L001_R1_001.good.fq -out=output/function/Konzo36_S74_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo36_S74_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo37_S73_ME_L001_R1_001.good.fq -out=output/function/Konzo37_S73_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo37_S73_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo38_S72_ME_L001_R1_001.good.fq -out=output/function/Konzo38_S72_ME_L001_R1.human_trimmed.fq -human=output/kraken2_output/Konzo38_S72_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R1_001.good.fq -out=output/function/L0049-33-10-76-Vilain-180822_S33_R1.human_trimmed.fq -human=output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R1_001.good.fq -out=output/function/L0060-44-22-38-Vilain-180822_S44_R1.human_trimmed.fq -human=output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R1_001.good.fq -out=output/function/L0063-47-10-81-Vilain-180822_S47_R1.human_trimmed.fq -human=output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R1_001.good.fq -out=output/function/L0067-51-10-76-Vilain-180822_S51_R1.human_trimmed.fq -human=output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L337-A34-Vilain-181127_S32_R1_001.good.fq -out=output/function/L337-A34-Vilain-181127_S32_R1.human_trimmed.fq -human=output/kraken2_output/L337-A34-Vilain-181127_S32.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L338-A20-Vilain-181127_S33_R1_001.good.fq -out=output/function/L338-A20-Vilain-181127_S33_R1.human_trimmed.fq -human=output/kraken2_output/L338-A20-Vilain-181127_S33.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L339-A31-Vilain-181127_S34_R1_001.good.fq -out=output/function/L339-A31-Vilain-181127_S34_R1.human_trimmed.fq -human=output/kraken2_output/L339-A31-Vilain-181127_S34.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L340-A33-Vilain-181127_S35_R1_001.good.fq -out=output/function/L340-A33-Vilain-181127_S35_R1.human_trimmed.fq -human=output/kraken2_output/L340-A33-Vilain-181127_S35.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L341-A26-Vilain-181127_S36_R1_001.good.fq -out=output/function/L341-A26-Vilain-181127_S36_R1.human_trimmed.fq -human=output/kraken2_output/L341-A26-Vilain-181127_S36.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L342-A13-Vilain-181127_S37_R1_001.good.fq -out=output/function/L342-A13-Vilain-181127_S37_R1.human_trimmed.fq -human=output/kraken2_output/L342-A13-Vilain-181127_S37.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L343-A32-Vilain-181127_S38_R1_001.good.fq -out=output/function/L343-A32-Vilain-181127_S38_R1.human_trimmed.fq -human=output/kraken2_output/L343-A32-Vilain-181127_S38.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L344-A3-Vilain-181127_S39_R1_001.good.fq -out=output/function/L344-A3-Vilain-181127_S39_R1.human_trimmed.fq -human=output/kraken2_output/L344-A3-Vilain-181127_S39.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L345-A6-Vilain-181127_S40_R1_001.good.fq -out=output/function/L345-A6-Vilain-181127_S40_R1.human_trimmed.fq -human=output/kraken2_output/L345-A6-Vilain-181127_S40.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0050-34-22-38-Vilain-180822_S34_R2_001.good.fq -out=output/function/L0050-34-22-38-Vilain-180822_S34_R2.human_trimmed.fq -human=output/kraken2_output/L0050-34-22-38-Vilain-180822_S34.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0070-54-10-079-Vilain-180822_S54_R2_001.good.fq -out=output/function/L0070-54-10-079-Vilain-180822_S54_R2.human_trimmed.fq -human=output/kraken2_output/L0070-54-10-079-Vilain-180822_S54.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0081-65-10-079-Vilain-180822_S65_R2_001.good.fq -out=output/function/L0081-65-10-079-Vilain-180822_S65_R2.human_trimmed.fq -human=output/kraken2_output/L0081-65-10-079-Vilain-180822_S65.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo10_S91_ME_L001_R2_001.good.fq -out=output/function/Konzo10_S91_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo10_S91_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo11_S90_ME_L001_R2_001.good.fq -out=output/function/Konzo11_S90_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo11_S90_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo12_S89_ME_L001_R2_001.good.fq -out=output/function/Konzo12_S89_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo12_S89_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo13_S88_ME_L001_R2_001.good.fq -out=output/function/Konzo13_S88_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo13_S88_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo14_S2_ME_L001_R2_001.good.fq -out=output/function/Konzo14_S2_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo14_S2_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo15_S3_ME_L001_R2_001.good.fq -out=output/function/Konzo15_S3_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo15_S3_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo16_S4_ME_L001_R2_001.good.fq -out=output/function/Konzo16_S4_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo16_S4_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo17_S5_ME_L001_R2_001.good.fq -out=output/function/Konzo17_S5_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo17_S5_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo18_S87_ME_L001_R2_001.good.fq -out=output/function/Konzo18_S87_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo18_S87_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo19_S86_ME_L001_R2_001.good.fq -out=output/function/Konzo19_S86_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo19_S86_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo1_S1_ME_L001_R2_001.good.fq -out=output/function/Konzo1_S1_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo1_S1_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo20_S85_ME_L001_R2_001.good.fq -out=output/function/Konzo20_S85_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo20_S85_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo21_S84_ME_L001_R2_001.good.fq -out=output/function/Konzo21_S84_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo21_S84_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo22_S83_ME_L001_R2_001.good.fq -out=output/function/Konzo22_S83_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo22_S83_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo23_S82_ME_L001_R2_001.good.fq -out=output/function/Konzo23_S82_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo23_S82_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo24_S81_ME_L001_R2_001.good.fq -out=output/function/Konzo24_S81_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo24_S81_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo25_S80_ME_L001_R2_001.good.fq -out=output/function/Konzo25_S80_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo25_S80_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo26_S7_ME_L001_R2_001.good.fq -out=output/function/Konzo26_S7_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo26_S7_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo27_S8_ME_L001_R2_001.good.fq -out=output/function/Konzo27_S8_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo27_S8_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo28_S9_ME_L001_R2_001.good.fq -out=output/function/Konzo28_S9_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo28_S9_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo29_S10_ME_L001_R2_001.good.fq -out=output/function/Konzo29_S10_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo29_S10_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo2_S6_ME_L001_R2_001.good.fq -out=output/function/Konzo2_S6_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo2_S6_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo30_S79_ME_L001_R2_001.good.fq -out=output/function/Konzo30_S79_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo30_S79_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo31_S78_ME_L001_R2_001.good.fq -out=output/function/Konzo31_S78_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo31_S78_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo32_S12_ME_L001_R2_001.good.fq -out=output/function/Konzo32_S12_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo32_S12_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo33_S77_ME_L001_R2_001.good.fq -out=output/function/Konzo33_S77_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo33_S77_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo34_S76_ME_L001_R2_001.good.fq -out=output/function/Konzo34_S76_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo34_S76_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo35_S75_ME_L001_R2_001.good.fq -out=output/function/Konzo35_S75_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo35_S75_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo36_S74_ME_L001_R2_001.good.fq -out=output/function/Konzo36_S74_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo36_S74_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo37_S73_ME_L001_R2_001.good.fq -out=output/function/Konzo37_S73_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo37_S73_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/Konzo38_S72_ME_L001_R2_001.good.fq -out=output/function/Konzo38_S72_ME_L001_R2.human_trimmed.fq -human=output/kraken2_output/Konzo38_S72_ME_L001.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0049-33-10-76-Vilain-180822_S33_R2_001.good.fq -out=output/function/L0049-33-10-76-Vilain-180822_S33_R2.human_trimmed.fq -human=output/kraken2_output/L0049-33-10-76-Vilain-180822_S33.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0060-44-22-38-Vilain-180822_S44_R2_001.good.fq -out=output/function/L0060-44-22-38-Vilain-180822_S44_R2.human_trimmed.fq -human=output/kraken2_output/L0060-44-22-38-Vilain-180822_S44.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0063-47-10-81-Vilain-180822_S47_R2_001.good.fq -out=output/function/L0063-47-10-81-Vilain-180822_S47_R2.human_trimmed.fq -human=output/kraken2_output/L0063-47-10-81-Vilain-180822_S47.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L0067-51-10-76-Vilain-180822_S51_R2_001.good.fq -out=output/function/L0067-51-10-76-Vilain-180822_S51_R2.human_trimmed.fq -human=output/kraken2_output/L0067-51-10-76-Vilain-180822_S51.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L337-A34-Vilain-181127_S32_R2_001.good.fq -out=output/function/L337-A34-Vilain-181127_S32_R2.human_trimmed.fq -human=output/kraken2_output/L337-A34-Vilain-181127_S32.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L338-A20-Vilain-181127_S33_R2_001.good.fq -out=output/function/L338-A20-Vilain-181127_S33_R2.human_trimmed.fq -human=output/kraken2_output/L338-A20-Vilain-181127_S33.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L339-A31-Vilain-181127_S34_R2_001.good.fq -out=output/function/L339-A31-Vilain-181127_S34_R2.human_trimmed.fq -human=output/kraken2_output/L339-A31-Vilain-181127_S34.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L340-A33-Vilain-181127_S35_R2_001.good.fq -out=output/function/L340-A33-Vilain-181127_S35_R2.human_trimmed.fq -human=output/kraken2_output/L340-A33-Vilain-181127_S35.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L341-A26-Vilain-181127_S36_R2_001.good.fq -out=output/function/L341-A26-Vilain-181127_S36_R2.human_trimmed.fq -human=output/kraken2_output/L341-A26-Vilain-181127_S36.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L342-A13-Vilain-181127_S37_R2_001.good.fq -out=output/function/L342-A13-Vilain-181127_S37_R2.human_trimmed.fq -human=output/kraken2_output/L342-A13-Vilain-181127_S37.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L343-A32-Vilain-181127_S38_R2_001.good.fq -out=output/function/L343-A32-Vilain-181127_S38_R2.human_trimmed.fq -human=output/kraken2_output/L343-A32-Vilain-181127_S38.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L344-A3-Vilain-181127_S39_R2_001.good.fq -out=output/function/L344-A3-Vilain-181127_S39_R2.human_trimmed.fq -human=output/kraken2_output/L344-A3-Vilain-181127_S39.non_human_reads.txt;
perl scripts/filter_fastq.pl -R=preprocessing/good/L345-A6-Vilain-181127_S40_R2_001.good.fq -out=output/function/L345-A6-Vilain-181127_S40_R2.human_trimmed.fq -human=output/kraken2_output/L345-A6-Vilain-181127_S40.non_human_reads.txt;
```

# Alignement agaisnt KEGG microbial genes with bowtie2/2.3.5.1

```
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0050-34-22-38-Vilain-180822_S34_R1.human_trimmed.fq  -2 output/function/L0050-34-22-38-Vilain-180822_S34_R2.human_trimmed.fq -S output/function/mapping/L0050-34-22-38-Vilain-180822_S34.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0070-54-10-079-Vilain-180822_S54_R1.human_trimmed.fq  -2 output/function/L0070-54-10-079-Vilain-180822_S54_R2.human_trimmed.fq -S output/function/mapping/L0070-54-10-079-Vilain-180822_S54.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0081-65-10-079-Vilain-180822_S65_R1.human_trimmed.fq  -2 output/function/L0081-65-10-079-Vilain-180822_S65_R2.human_trimmed.fq -S output/function/mapping/L0081-65-10-079-Vilain-180822_S65.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo10_S91_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo10_S91_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo10_S91_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo11_S90_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo11_S90_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo11_S90_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo12_S89_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo12_S89_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo12_S89_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo13_S88_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo13_S88_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo13_S88_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo14_S2_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo14_S2_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo14_S2_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo15_S3_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo15_S3_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo15_S3_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo16_S4_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo16_S4_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo16_S4_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo17_S5_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo17_S5_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo17_S5_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo18_S87_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo18_S87_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo18_S87_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo19_S86_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo19_S86_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo19_S86_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo1_S1_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo1_S1_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo1_S1_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo20_S85_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo20_S85_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo20_S85_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo21_S84_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo21_S84_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo21_S84_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo22_S83_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo22_S83_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo22_S83_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo23_S82_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo23_S82_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo23_S82_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo24_S81_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo24_S81_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo24_S81_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo25_S80_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo25_S80_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo25_S80_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo26_S7_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo26_S7_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo26_S7_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo27_S8_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo27_S8_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo27_S8_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo28_S9_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo28_S9_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo28_S9_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo29_S10_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo29_S10_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo29_S10_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo2_S6_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo2_S6_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo2_S6_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo30_S79_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo30_S79_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo30_S79_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo31_S78_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo31_S78_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo31_S78_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo32_S12_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo32_S12_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo32_S12_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo33_S77_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo33_S77_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo33_S77_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo34_S76_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo34_S76_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo34_S76_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo35_S75_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo35_S75_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo35_S75_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo36_S74_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo36_S74_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo36_S74_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo37_S73_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo37_S73_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo37_S73_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/Konzo38_S72_ME_L001_R1.human_trimmed.fq  -2 output/function/Konzo38_S72_ME_L001_R2.human_trimmed.fq -S output/function/mapping/Konzo38_S72_ME_L001.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0049-33-10-76-Vilain-180822_S33_R1.human_trimmed.fq  -2 output/function/L0049-33-10-76-Vilain-180822_S33_R2.human_trimmed.fq -S output/function/mapping/L0049-33-10-76-Vilain-180822_S33.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0060-44-22-38-Vilain-180822_S44_R1.human_trimmed.fq  -2 output/function/L0060-44-22-38-Vilain-180822_S44_R2.human_trimmed.fq -S output/function/mapping/L0060-44-22-38-Vilain-180822_S44.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0063-47-10-81-Vilain-180822_S47_R1.human_trimmed.fq  -2 output/function/L0063-47-10-81-Vilain-180822_S47_R2.human_trimmed.fq -S output/function/mapping/L0063-47-10-81-Vilain-180822_S47.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L0067-51-10-76-Vilain-180822_S51_R1.human_trimmed.fq  -2 output/function/L0067-51-10-76-Vilain-180822_S51_R2.human_trimmed.fq -S output/function/mapping/L0067-51-10-76-Vilain-180822_S51.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L337-A34-Vilain-181127_S32_R1.human_trimmed.fq  -2 output/function/L337-A34-Vilain-181127_S32_R2.human_trimmed.fq -S output/function/mapping/L337-A34-Vilain-181127_S32.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L338-A20-Vilain-181127_S33_R1.human_trimmed.fq  -2 output/function/L338-A20-Vilain-181127_S33_R2.human_trimmed.fq -S output/function/mapping/L338-A20-Vilain-181127_S33.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L339-A31-Vilain-181127_S34_R1.human_trimmed.fq  -2 output/function/L339-A31-Vilain-181127_S34_R2.human_trimmed.fq -S output/function/mapping/L339-A31-Vilain-181127_S34.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L340-A33-Vilain-181127_S35_R1.human_trimmed.fq  -2 output/function/L340-A33-Vilain-181127_S35_R2.human_trimmed.fq -S output/function/mapping/L340-A33-Vilain-181127_S35.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L341-A26-Vilain-181127_S36_R1.human_trimmed.fq  -2 output/function/L341-A26-Vilain-181127_S36_R2.human_trimmed.fq -S output/function/mapping/L341-A26-Vilain-181127_S36.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L342-A13-Vilain-181127_S37_R1.human_trimmed.fq  -2 output/function/L342-A13-Vilain-181127_S37_R2.human_trimmed.fq -S output/function/mapping/L342-A13-Vilain-181127_S37.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L343-A32-Vilain-181127_S38_R1.human_trimmed.fq  -2 output/function/L343-A32-Vilain-181127_S38_R2.human_trimmed.fq -S output/function/mapping/L343-A32-Vilain-181127_S38.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L344-A3-Vilain-181127_S39_R1.human_trimmed.fq  -2 output/function/L344-A3-Vilain-181127_S39_R2.human_trimmed.fq -S output/function/mapping/L344-A3-Vilain-181127_S39.kegg_mapping.sam
bowtie2 --no-unal --no-sq --omit-sec-seq --threads 32 -x db/bowtie2_keggDB -1 output/function/L345-A6-Vilain-181127_S40_R1.human_trimmed.fq  -2 output/function/L345-A6-Vilain-181127_S40_R2.human_trimmed.fq -S output/function/mapping/L345-A6-Vilain-181127_S40.kegg_mapping.sam
```

# Parse SAM alignement results

```
perl scripts/parse_sam.pl -query=output/function/mapping/L0050-34-22-38-Vilain-180822_S34.kegg_mapping.sam -out=output/function/mapping/L0050-34-22-38-Vilain-180822_S34.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0070-54-10-079-Vilain-180822_S54.kegg_mapping.sam -out=output/function/mapping/L0070-54-10-079-Vilain-180822_S54.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0081-65-10-079-Vilain-180822_S65.kegg_mapping.sam -out=output/function/mapping/L0081-65-10-079-Vilain-180822_S65.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo10_S91_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo10_S91_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo11_S90_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo11_S90_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo12_S89_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo12_S89_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo13_S88_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo13_S88_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo14_S2_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo14_S2_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo15_S3_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo15_S3_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo16_S4_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo16_S4_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo17_S5_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo17_S5_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo18_S87_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo18_S87_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo19_S86_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo19_S86_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo1_S1_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo1_S1_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo20_S85_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo20_S85_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo21_S84_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo21_S84_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo22_S83_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo22_S83_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo23_S82_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo23_S82_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo24_S81_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo24_S81_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo25_S80_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo25_S80_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo26_S7_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo26_S7_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo27_S8_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo27_S8_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo28_S9_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo28_S9_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo29_S10_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo29_S10_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo2_S6_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo2_S6_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo30_S79_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo30_S79_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo31_S78_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo31_S78_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo32_S12_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo32_S12_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo33_S77_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo33_S77_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo34_S76_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo34_S76_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo35_S75_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo35_S75_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo36_S74_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo36_S74_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo37_S73_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo37_S73_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/Konzo38_S72_ME_L001.kegg_mapping.sam -out=output/function/mapping/Konzo38_S72_ME_L001.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0049-33-10-76-Vilain-180822_S33.kegg_mapping.sam -out=output/function/mapping/L0049-33-10-76-Vilain-180822_S33.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0060-44-22-38-Vilain-180822_S44.kegg_mapping.sam -out=output/function/mapping/L0060-44-22-38-Vilain-180822_S44.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0063-47-10-81-Vilain-180822_S47.kegg_mapping.sam -out=output/function/mapping/L0063-47-10-81-Vilain-180822_S47.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L0067-51-10-76-Vilain-180822_S51.kegg_mapping.sam -out=output/function/mapping/L0067-51-10-76-Vilain-180822_S51.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L337-A34-Vilain-181127_S32.kegg_mapping.sam -out=output/function/mapping/L337-A34-Vilain-181127_S32.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L338-A20-Vilain-181127_S33.kegg_mapping.sam -out=output/function/mapping/L338-A20-Vilain-181127_S33.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L339-A31-Vilain-181127_S34.kegg_mapping.sam -out=output/function/mapping/L339-A31-Vilain-181127_S34.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L340-A33-Vilain-181127_S35.kegg_mapping.sam -out=output/function/mapping/L340-A33-Vilain-181127_S35.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L341-A26-Vilain-181127_S36.kegg_mapping.sam -out=output/function/mapping/L341-A26-Vilain-181127_S36.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L342-A13-Vilain-181127_S37.kegg_mapping.sam -out=output/function/mapping/L342-A13-Vilain-181127_S37.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L343-A32-Vilain-181127_S38.kegg_mapping.sam -out=output/function/mapping/L343-A32-Vilain-181127_S38.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L344-A3-Vilain-181127_S39.kegg_mapping.sam -out=output/function/mapping/L344-A3-Vilain-181127_S39.kegg_mapping.filtered.sam;
perl scripts/parse_sam.pl -query=output/function/mapping/L345-A6-Vilain-181127_S40.kegg_mapping.sam -out=output/function/mapping/L345-A6-Vilain-181127_S40.kegg_mapping.filtered.sam;
```

# concatenation of alignment into KEGG ko, module , pathway matrices

```
perl scripts/functional_annotation_metagenomic.pl -query=output/function/list.txt -ko=scripts/ko_genes.list -module=scripts/ko_module.list -pathway=scripts/ko_pathway.list -out_ko=output/function/konzo.dna_kegg_metaG.ko.count.txt -out_mod=output/function/konzo.dna_kegg_metaG.module.count.txt -out_path=output/function/konzo.dna_kegg_metaG.pathway.count.txt
```
