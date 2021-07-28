#!/usr/bin/perl -w
use strict;
use Getopt::Long;
###############################
# CONCATENATION OF centrifuge_report
# perl concat.pl -query=[listof sample files] -output=[normalized data output] -out_spe=[normalized species data output] -out_raw=[unique raw count data output] -out_rawsp=[unique raw count species data output]
###############################
my $usage = "Erreur!!  Usage: perl concat.pl -query=[listof sample files] -output=[normalized data output] -out_spe=[normalized species data output] -out_raw=[unique raw count data output] -out_rawsp=[unique raw count species data output] $!\n";
my $options = join(" ", @ARGV);
my ($query, $output, $output2,$output3, $map, $kof, $mof);
GetOptions ("query=s"               => \$query,
            "pathway=s"             => \$map,
            "ko=s"                  => \$kof,
            "module=s"                  => \$mof,
            "out_ko=s"              => \$output,
            "out_path=s"              => \$output2,
            "out_mod=s"              => \$output3,
            ) or die("Error in command line arguments\n");
###############################################################################
my (%hash, %gene2ko, %ko2path,%ko2mod, %ko_out, %path_out, %mod_out, %mod_out2);
###############################################################################
open(FILE, $kof) or die "can't open $kof";
while (<FILE>) {
    chomp;
    my @tab = split ("\t", $_);
    my $dd;
    if ($tab[0] =~ /ko:(.*)$/){
      $dd = $1;
    }
    $gene2ko{$tab[1]}{$1} = 1;
}
close FILE;
###############################################################################
my (%mod, %modcount);
open(FILE, $mof) or die "can't open $mof";
while (<FILE>) {
  chomp;
  my @tab = split ("\t", $_);
  if ($tab[1] =~ /^md:/){
		if ($tab[0] =~ /ko:(.*)$/){
			$ko2mod{$1}{$tab[1]} = 1;
      if (! $modcount{$tab[1]}{$1}){
        $modcount{$tab[1]}{$1} =1;
        $mod{$tab[1]}++;
      }
		}
	}
}
close FILE;
###############
###############################################################################
open(FILE, $map) or die "can't open $map";
while (<FILE>) {
  chomp;
  my @tab = split ("\t", $_);
  if ($tab[1] =~ /^path:ko/){
		if ($tab[0] =~ /ko:(.*)$/){
			$ko2path{$1}{$tab[1]} = 1;
		}
	}
}
close FILE;
###############################################################################
my @sample;
###############################################################################
open(FILE, $query) or die "can't open $query";
while (<FILE>) {
    chomp;
    my $s;
    if (($_ =~ /.*\/(.*)\.kegg_mapping\.filtered\.sam/)){
      $s = $1;
      push (@sample, $s);
    }
    else{
      print "don't have the right regex for file list parsing: $_.\nstopped process\nSearch for *.kegg_mapping\.sam files, Change file name or script parser\n";
      exit;
    }
    open(FILEX, $_) or die "can't open $_";
    my %read;
		while (<FILEX>) {
      chomp;
      my @tab = split ("\t", $_);
      if ($gene2ko{$tab[1]}){
        foreach my $i (keys %{$gene2ko{$tab[1]}}){
            if (! $read{$tab[0]}){
							if (! $ko_out{$i}{$s}) {
								$ko_out{$i}{$s} = 1 ;
							  $read{$tab[0]} = 1;
								if ($ko2path{$i}){
									foreach my $i2 (keys %{$ko2path{$i}})	{
										#print $i."\t".$i2."\n";
										#my $u = <stdin>;
										if (! $path_out{$i2}{$s}){
									     $path_out{$i2}{$s} = 1 ;
									  }
									  else{
									     $path_out{$i2}{$s}++;
									  }
									}
								}
                if ($ko2mod{$i}){
									foreach my $i2 (keys %{$ko2mod{$i}})	{
										#print $i."\t".$i2."\n";
										#my $u = <stdin>;
										if (! $mod_out{$i2}{$s}) {
									     $mod_out{$i2}{$s} = 1 ;
                       $mod_out2{$i2}{$s}{$i} =1;
									  }
										else{
									     $mod_out{$i2}{$s}++;
                       $mod_out2{$i2}{$s}{$i} =1;
									  }
									}
								}
						}
						else{
								$ko_out{$i}{$s}++;
								$read{$tab[0]} = 1;
								if ($ko2path{$i}){
                   foreach my $i2 (keys %{$ko2path{$i}}) {
                     if (! $path_out{$i2}{$s}){
                        $path_out{$i2}{$s} = 1 ;
                     }
                     else{
                        $path_out{$i2}{$s}++;
                     }
										 #print $i."\t".$i2."\n";
										 #my $u = <stdin>;
                   }
                 }
                 if ($ko2mod{$i}){
 									foreach my $i2 (keys %{$ko2mod{$i}})	{
 										#print $i."\t".$i2."\n";
 										#my $u = <stdin>;
 										if ((! $mod_out{$i2}{$s}) && (! $mod_out2{$i2}{$s}{$i})){
 									     $mod_out{$i2}{$s} = 1 ;
                        $mod_out2{$i2}{$s}{$i} =1;
 									  }
 									  elsif (! $mod_out2{$i2}{$s}{$i}){
 									     $mod_out{$i2}{$s}++;
                        $mod_out2{$i2}{$s}{$i} =1;
 									  }
 									}
 								}
							}
					 }
        }
      }
    }
    close FILEX;
    close W;
}
close FILE;

#foreach my $i (keys %path_out){
#	foreach my $i2 (keys %{$path_out{$i}}){
#		print $i."\t".$i2."\n";
#		my $u = <stdin>;
#	}
#}


open(W, ">$output") or die "pbm output $usage";
my $sa = join("\t", @sample);
print W "ko\t$sa\n";
foreach my $i (sort keys %ko_out){
  print W $i;
  for my $i2 (0..$#sample){
    if ($ko_out{$i}{$sample[$i2]}){
      print W "\t".$ko_out{$i}{$sample[$i2]};
    }
    else{
      print W "\t0";
    }
  }
  print W "\n";
}
close W;
open(W2, ">$output2") or die "pbm output $usage";
print W2 "pathway\t$sa\n";
foreach my $i (sort keys %path_out){
  print W2 $i;
  for my $i2 (0..$#sample){
    if ($path_out{$i}{$sample[$i2]}){
      print W2 "\t".$path_out{$i}{$sample[$i2]};
    }
    else{
      print W2 "\t0";
    }
  }
  print W2 "\n";
}
close W2;
open(W3, ">$output3") or die "pbm output $usage";
print W3 "module\t$sa\n";
foreach my $i (sort keys %mod_out){
  print W3 $i;
  for my $i2 (0..$#sample){
    if ($mod_out{$i}{$sample[$i2]}){
      print W3 "\t".$mod_out{$i}{$sample[$i2]};
    }
    else{
      print W3 "\t0";
    }
  }
  print W3 "\n";
}
close W3;
