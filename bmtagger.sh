#!/bin/bash
export PATH=/home/nvashist/yes/bin/:$PATH

for f in /data/Neerja/Konzo3_MetagenomicMicrobiome/raw_fastq/20181127_Vilain_Konzo3_part1of2/FASTQ_Generation_2019-01-26_03_17_35Z-153995021/*/
do bmtagger.sh -C /data/Neerja/scripts/bmtagger.conf -b /data/Neerja/hg38/Homo_sapiens_assembly38.bitmask -x /data/Neerja/hg38/Homo_sapiens_assembly38.srprism -T /data/Neerja/temp -q 1 -1 $f/*R1*.fastq -2 $f/*R2*.fastq -X -o /data/Neerja/Konzo3_output/bmtagger/$(basename $f)
done
for f in /data/Neerja/Konzo3_MetagenomicMicrobiome/raw_fastq/Konzo3_Part2of2/Vilain_Konzo3_2of2-116664549/FASTQ_Generation_2019-02-09_02_08_07Z-159815670/*/
do bmtagger.sh -C /data/Neerja/scripts/bmtagger.conf -b /data/Neerja/hg38/Homo_sapiens_assembly38.bitmask -x /data/Neerja/hg38/Homo_sapiens_assembly38.srprism -T /data/Neerja/temp -q 1 -1 $f/*R1*.fastq -2 $f/*R2*.fastq -X -o /data/Neerja/Konzo3_output/bmtagger/$(basename $f)
done
