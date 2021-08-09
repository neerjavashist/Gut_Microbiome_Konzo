#!/usr/bin/perl -w
#
use strict;
use Getopt::Long;

my $usage = "Erreur!!  Usage: $0 $!\n";

###############################################################################

# PARAMETERS

my $NA= "NA";
my ($queryfile,$queryfile2, $out,$out2, $theo);
my $options = join(" ", @ARGV);

GetOptions ("R=s"               => \$queryfile,
            "out=s"              => \$out,
            "human=s"              => \$theo,
            ) or die("Error in command line arguments\n");
###############################################################################
###############################################################################
print localtime()."\t"."charging non human readIDs \n";
my (%hash);
open(F, $theo) or die " No non human file $usage";
while (<F>) {
    chomp;
    my @tab = split("\t", $_);
    $hash{$tab[0]} = 1;
}
close F;

###############################################################################
print localtime()."\t"."parsing and writing $queryfile\n";
my $astuce = 2;
open(W, ">$out") or die "can't open $out\n";
open(F, $queryfile) or die "No  $queryfile $usage";
while (<F>) {
    chomp;
    if ($_ =~ /^@(\S+)/){
      if ($hash{$1}){
        $astuce = 4;
        print W $_."\n";
      }
      else{
        $astuce = 2;
      }
    }
    elsif ($astuce eq "4"){
      print W $_."\n";
    }
}
close F;
close W;
###############################################################################
