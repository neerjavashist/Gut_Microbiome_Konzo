#!/bin/bash
export PATH=/home/nvashist/Bracken/:$PATH
INDIR=/data/Neerja/Konzo3_output/kraken2
OUTDIR=/data/Neerja/Konzo3_output/Bracken
for f in {256..407}
do
bracken -d /data/Neerja/Kraken_DB_0303 -i $INDIR/L$f.Kraken_report -o $OUTDIR/Konzo_$f.bracken_report -r 150
done
