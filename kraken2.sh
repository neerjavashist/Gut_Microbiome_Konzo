#!/bin/bash
export PATH=/home/nvashist/kraken2/:$PATH

INDIR=/data/Neerja/Konzo3_output/skewer
OUTDIR=/data/Neerja/Konzo3_output/kraken2
for f in {256..407}
do
Fastq1=`ls $INDIR/L$f*_*.fastq.gz | tr "\n" " " `

kraken2 --db /data/Neerja/Kraken_DB_0303 --paired --use-names \
--unclassified-out $OUTDIR/reads/L$f.unclassifed# \
--classified-out $OUTDIR/reads/L$f.classifed# \
--output - \
--report-zero-counts \
--report $OUTDIR/L$f.Kraken_report \
$Fastq1 2> $OUTDIR/L$f.Kraken_out


kraken2 --db /data/Neerja/Kraken_DB_0303 --paired --use-names \
--unclassified-out $OUTDIR/reads/L$f.unclassifed# \
--classified-out $OUTDIR/reads/L$f.classifed# \
--output - \
--use-mpa-style \
--report $OUTDIR/L$f.Kraken_mpa \
$Fastq1 &&
awk '{print$0}' $OUTDIR/L$f.Kraken_mpa > $OUTDIR/L$f.metagenomics.txt
done
