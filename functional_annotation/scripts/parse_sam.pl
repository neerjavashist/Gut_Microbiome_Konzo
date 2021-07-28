#!/usr/bin/perl -w
#
use strict;
use Getopt::Long;

my $usage = "Erreur!!  Usage: $0 $!\n";

###############################################################################

# CATALOGUE DATA
my $lenfile;

# NCBI / REFSEQ
my $refseq;
my $nodefile;
my $namefile;
my $mergedfile;
# PARAMETERS

my $NA= "NA";

my ($queryfile, $blastkegg, $blastrefseq, $out,$out2,$out3,$out4,$out5,$out6, $theo);
my $options = join(" ", @ARGV);

GetOptions ("query=s"               => \$queryfile,
            "out=s"              => \$out,
            ) or die("Error in command line arguments\n");
###############################################################################
###############################################################################
###############################################################################
# RETRIEVE NCBI TAXONOMY
###############################################################################

my (%hash,%hash2, %mapq, $ok, %count);
open (W, ">$out") or die "can't open $out\n";
open(F, $queryfile) or die "No Taxonomy queryfile $usage";
while (<F>) {
    chomp;
    if ($_ !~ /^@/){
      my @tab = split("\t", $_);
      if (! $hash{$tab[0]}){
        $hash{$tab[0]} = $tab[2];
        $mapq{$tab[0]} = $tab[4];
      }
      elsif ($hash{$tab[0]} eq $tab[2]){
        $ok++;
      }
      elsif ($tab[4] > $mapq{$tab[0]} ){
        $hash2{$tab[0]}{$tab[2]} = $tab[1];
        $hash{$tab[0]} = $tab[2];
      }
      else{
        $hash2{$tab[0]}{$tab[2]} = $tab[1];
      }
    }
}
close F;

foreach my $i (sort keys %hash){
  print W $i."\t".$hash{$i}."\n";
}

print "Number of accorded pairs: $ok\n";
