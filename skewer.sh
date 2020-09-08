#!/bin/bash
export PATH=/home/nvashist/yes/bin/:$PATH
for f in /data/Neerja/Konzo3_output/bmtagger/*_1.fastq
do skewer -l 30 -m pe -x CTGTCTCTTATACACATCT -y CTGTCTCTTATACACATCT -z -o /data/Neerja/Konzo3_output/skewer/$(basename $f .fastq) $f /data/Neerja/Konzo3_output/bmtagger/$(basename $f _1.fastq)*_2.fastq
done
