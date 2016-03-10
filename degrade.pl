#!/usr/bin/perl
#
# Last edit: 7-Feb-15
# Author: Tyler K. Chafin - tkchafin@uark.edu
# Title: fragment.pl 
# This script was developed in order to simulate random shearing 
# of input DNA sequences. This is accomplished by randomly sampling 
# target fragment lengths from a distribution of user-defined parameters. 
# User may choose to sample these fragment lengths from either a uniform
# or Gaussian distribution, bounded or unbounded. 
#
#FIX : Program hangs at small mean; not sure whats up
#
use strict; 
use warnings; 
use Getopt::Long;

#DEBUG
my $mem=0;


#Initialize variables
my @input; 
my $dist = "g"; 
my @frags; 
my $index = 0; 
my $first = 0; 
my @rand; 
my $mean = 1000;
my $stdev = 1000; 
my $l_bound = 0; 
my $u_bound = 1000000000; 
my $samples= 1000000;
my $before = 0; 
my @newfrags;
my $outseqs=0; 
my $help = 0; 
my $outname = "fragment";

#Call command-line parsing subroutine
parseArgs();


#Some checks for valid option values...
$dist =~ /^[^gu]/i and die "\nPROGRAM KILLED:\nOh no! You seem to have selected a non-existant distribution option! Try again using the \"-d\" option. \n\n"; 
if ($dist =~ /^n/i){ 
	$l_bound > $mean and print "\nWARNING:\nYou have set the lower bound for your sampled fragment size distribution to higher than the mean!\n\n";
	$u_bound < $mean and print "\nWARNING:\nYou have set the upper bound for your sampled fragment size distribution to lower than the mean!\n\n";
}

########################################################
#Sample selected distribution type to create target fragment sizes

#Check if distribution type chosen is normal
if ($dist =~ /^g/i){
	#Sample from a bounded Gaussian distribution
	@rand = @{bounded_gaussian_rand($samples, $mean, $stdev, $l_bound, $u_bound)};
}

#Make another for a uniform distribution
if ($dist =~ /^u/i){ 
	@rand = @{bounded_uniform_rand($samples, $l_bound, $u_bound)};
}

########################################################

#Capture all "genomic" fragments (e.g. provided contigs/ scaffolded contigs)

@input = glob "@input";
#print "@input\n";

#Loop through each file; create array of all fragments
foreach my $file (@input){ 
	open (FILE, "$file") || die "Cannot open file $file: $!\n";
		while (<FILE>){ 
			chomp; 
			if ($_ =~ /^\>/){ 	
				#print "$_\n";		
				if ($first == 0){ 
					$first++;
					next;
				}else{ 
					$index++; 
					next;
				} 
			}else{ 
				#print "$_\n";
				if ($frags[$index]){
					$frags[$index] .= $_;
				}else{ 
					$frags[$index] = ""; 
					$frags[$index] .= $_;
				}
			}
		}
}

if ($before == 1){ 
	#Call subroutine fragtable()
	fragtable(\@frags, "pre_fragment"); 
}

		
####################################################


#Random fragmentation according to selected size distribution (e.g. uniform or normal)

#Initial values for indices and substr parameters
my $rind = 0; 
my $offset = 0; 
my $length = int($rand[$rind]); 
my $total = 0; 
my $full_length; 

foreach (@frags){ 
	
	#Foreach new fragment, reset total substr bases
	$total = 0; 
	$full_length = length($_); 
		print "Unfragmented length is $full_length\n";
	

#
#INFINITE LOOP UNDER CERTAIN CONDITIONS: NEEDS REVISION!!!!!!!
#
#
	
	#Keep sequentially substr fragment until it has all been fragmented
	while ($full_length >= $total){
	 	
		#If we have used all currently captured random samples, resample distribution of frag lengths
		if ($rind >= $samples){ 
			print "Resetting random number distribution...\n"; 
			$dist =~ /^g/i and @rand = @{bounded_gaussian_rand($samples, $mean, $stdev, $l_bound, $u_bound)};
			$dist =~ /^u/i and @rand = @{bounded_uniform_rand($samples, $l_bound, $u_bound)};
			#Reset values
			$rind = 0; 
			#Reset values
			$rind = 0; 
			$offset = 0;
			$length = int($rand[$rind]); 
		}
		
		if ($length >= $full_length){ 
			push (@newfrags, $_);
			last; 
		}
	#	
		$total += $length and print "$total\n"; 
			#print "Fragmented length is $total\n"; #Testing purposes
		$total > $full_length and $length -= ($total - $full_length); 
		$length <= 0 and next; 
		$offset > $full_length and last; 
		
		#Subtring fragment, capture sub in nefrags array
		my $sub = substr($_, $offset, $length); 
			#print $sub . "\n"; #Uncomment for testing purposes
		push (@newfrags, $sub) unless $mem==1;
 	#my $temp = length($sub);
		#print "$temp\n"; 
		#Set offset to current random length (e.g. start after end of last one)
		$offset = $total; 
		$length = int($rand[$rind+1]);
		#Increment index for random number array
		$rind++;  
	}
	
}

if ($mem==0){
print "Writing fragtable...\n";
fragtable(\@newfrags, $outname);
#print join("\n", @rand);

$outseqs == 1 and print "Writing FASTA file...\n" and fragFASTA(\@newfrags, $outname); 
}

exit;


#################################################################
######################## SUBROUTINES ############################
#################################################################


sub bounded_gaussian_rand {

#Modified from: 
#O'Rielly Perl Cookbook 
#By Tom Christiansen & Nathan Torkington; ISBN 1-56592-243-3, 794 pages.
#First Edition, August 1998. 
#Chapter 2.10 "Generating Biased Random Numbers"	

#Implements polar Box Muller method to pseudosample a Gaussian distribution

print "Bounded gaussian\n"; 
	my @rsamp; 
	my ($num, $m, $sd, $min, $max) = @_; 
	$num || die "bounded_gaussian_rand: Number of samples not specified!\n";
	$m || die "bounded_gaussian_rand: Mean for sampled distribution not specified!\n";
	$sd || die "bounded_gaussian_rand: Standard deviation not specified!\n";
	#You can choose to not specify bounds if you prefer
	
	#print "$min , $max\n\n";
	
	until(scalar(@rsamp) > $num){;
		my ($u1, $u2);  # uniformly distributed random numbers
		my $w;          # variance, then a weight
		my ($g1, $g2);  # gaussian-distributed numbers

		do {
			$u1 = 2 * rand() - 1;
			$u2 = 2 * rand() - 1;
			$w = $u1*$u1 + $u2*$u2;
		} while ( $w >= 1 );

		$w = sqrt( (-2 * log($w))  / $w );
		$g2 = $u1 * $w;
		$g1 = $u2 * $w;
		
		#Modify to fit new distribution
		$g1 = $g1 * $sd + $m;
		$g2 = $g2 * $sd + $m;
 

		#Check if number within bounds
		if ($min ne ""){ 
			$g1 >= $min and push (@rsamp, $g1);
			$g2 >= $min and push (@rsamp, $g2);
			next;
		}else{ 
			$max || push (@rsamp, $g1);
			$max || push (@rsamp, $g2);
			next;
		}
		if ($max){ 
			$g1 <= $max and push (@rsamp, $g1);
			$g2 <= $max and push (@rsamp, $g2);
			next; 
		}else{ 
			push @rsamp, $g1;
			push @rsamp, $g2;
			next;
		}
	}
	#print "@rsamp\n";
	return \@rsamp;
} 

########################################################################

sub bounded_uniform_rand{
		#Randomly sample uniform distribution and return.. Same format as
		# used above for Gaussian
		my @rsamp; 
		my  ($num, $min, $max) = @_; 
		my $r; 
		
		#Until @rsamp contains as many sampled numbers as we want, 
		# continue loop
		until (scalar(@rsamp) > $num){ 
			
			$r = rand($max); 
			$r < $min and next; 
			$r > $max and next; 
			push (@rsamp, $r); 			
		}
		return \@rsamp;
}

########################################################################

sub fragtable{ 
	#Recieve an array of strings and return a hash of frequency of each
	# string length (i.e. Key= length; value= # occurances). Return 
	my ($array, $out) = @_; 
	my %hash; 
	
	foreach(@{$array}){ 
		$hash{length($_)}++; 
	}
	open ( OUT, ">$out.tsv" ) || die "Can't open $out.tsv!\n\n";

	#Sort the lengths and print them, followed by the number of times they occured

    foreach (sort {$a <=> $b} keys %hash){
        print OUT "\n$_\t$hash{$_}";
    }
    
	close OUT;
}

########################################################################
#Simple subroutine to output fragment array to a single fasta file, 
# where the header is the fragment length
sub fragFASTA{ 
	my ($array, $out) = @_; 
	
	open (OUT, ">$out.fas") || die "Can't open $out.fas!\n\n"; 
	
	foreach (@{$array}){ 
		print OUT ">" . length($_) . "\n$_\n";
	}
	close OUT; 
}

########################################################################

sub parseArgs{ 
	
	my $usage = "\n\n####################
fragment.pl v.0.1

Author: Tyler K. Chafin - tkchafin\@uark.edu

Last Modified: 7 February 2015

This script was developed in order to simulate random shearing of input DNA sequences. This is accomplished by randomly sampling target fragment lengths from a distribution of user-defined parameters. User may choose to sample these fragment lengths from either a uniform or Gaussian distribution, bounded or unbounded. 
####################

Usage: $0 -i /path/to/*.fasta

Input should be in FASTA format. Headers will be ignored. Sequences will be digested as is, and not concatenated. 

Mandatory Arguments
	-i, --input	-  Path to input FASTA file(s) containing your genomic sequences [e.g. /home/*.fasta]
		
Optional Arguments
	-d, --dist	-  Distribution type, normal or uniform [default = normal]
	-o, --out	-  Prefix name for output [default = fragment]
	-u, --up	-  Upper bound for fragment size dist [default = none]
	-l, --low	-  Lower bound for fragment size dsit [default = none]
	-b, --pre	-  Toggle on (boolean) to print pre=fragmentation table 
	-f, --fas	-  Toggle on (boolean) to print final fragments to FASTA file
		
Mandatory if distribution type = normal 
	-m, --mean	-  Mean for fragment size distribution
	-s, --sd	-  Standard deviation for fragment size distribution
		
fragment.pl will by default output a tsv file containing fragment lengths for randomy sheared sequences, with one column for the fragment length and a second for the number of occurances. You can optionally also output sheared sequences to a FASTA file, as well as another tsv file for the unsheared DNA. 

";

	my $result = GetOptions
	( 
	'help|h!'		=> \$help, 
	'input|i=s{1,}'	=> \@input,
	'out|o=s'		=> \$outname,
	'b|pre!'		=> \$before, 
	'f|fas!'		=> \$outseqs, 
	'u|up=i'		=> \$u_bound, 
	'l|low=i'		=> \$l_bound, 
	'm|mean=i'		=> \$mean, 
	's|sd=i'		=> \$stdev, 
	'd|dist=s'		=> \$dist, 
	); 

$help == 1 and die "$usage\n";
@input or die "\nPROGRAM KILLED: Input not specified!\n$usage\n";
	
	
}
