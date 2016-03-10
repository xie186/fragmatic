#!/usr/bin/perl 

use strict; 
#use warnings; 
use Getopt::Long; 

our $infile;
our $up;
our $low; 


my $usage="\nUsage: $0 -i infile.tsv -u 330 -l 270

This script is intended to parse the output from ddRADsim.pl

	-i, --infile	- .tsv file output from ddRADsim.pl
	-u, --up 	- Upper limits for size selection range
	-l, --lower 	- Lower limit for size selection range

"; 

GetOptions(
'infile|i=s'	=> \$infile,
'up|u=i'	=> \$up,
'low|l=i'	=> \$low,
); 

$infile or die "\n$usage\n";
$up or die "\n$usage\n";
$low or die "\n$usage\n";

my @line; 
my $count=0; 

open (TSV, "$infile") || die "Can't open $infile!\n"; 

while (<TSV>){ 
    chomp; 
    @line = split /\t/, $_; 
    if ($line[0] <= $up && $line[0] >= $low){ 
	$count+= $line[1]; 
    }
}

close TSV; 

print "Number of fragments in size range $low  -  $up : $count\n"; 


exit; 
