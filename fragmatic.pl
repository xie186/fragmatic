#!/usr/bin/perl
#
# Last edit: 7-May-15
# Author: Tyler K. Chafin - tkchafin@uark.edu
# Algorithmic improvements by Steve M. Mussmann and Bradley T. Martin
# Title: fragsim.pl 
# This script was developed to provide a means to digest genomic DNA 
# sequences in silico.
#

use strict;
use warnings; 
use File::Basename;
use Getopt::Std;

my %opts;
getopts( 'i:r:o:fh', \%opts );

#parse command-line arguments
my ($input, $re, $outname, $fasta) = parseArgs( \%opts ); 

print "$re\n\n";

my @re = split /\s+/, $re;

#Initialize some local variables
my ( $filepath, $dirpath ) = fileparse ( $outname );
my @count;
my $len; 
my $key;  
my %total;
my $newstart;
my @frags; 
my @re_full;
my @re_pattern;
my @plainfreq; 
my %sites; 
my $count; 

#Hash for iupac codes
my %iupac =(

	"N" => ["A","T","G","C"],
	"R" => ["A", "G"],
	"Y" => ["T", "C"], 
	"K" => ["G", "T"], 
	"M" => ["A", "C"],
	"S" => ["G", "C"],
	"W" => ["A", "T"], 
	"B" => ["C", "G", "T"],
	"D" => ["A", "G", "T"],
	"H" => ["A", "C", "T"],
	"V" => ["A", "C", "G"],
	
);

########################################################################

#Add cut site identifiers to flank each fragment
#Capture cut site patterns for later
#Resolve degenerate bases in cut sites, if present

#Initialize some variables
my $re_num = 0;
my $re_full; 
my $re_pattern;  
my $temp; 
my @re_plain; 

#Convert RE sites to uppercase
$_ = uc for @re; 
#Maybe add some more printing of options called...
print "\nInput restriction sites: @re\n\n"; 
foreach my $site (@re){ 
	
	#Check that RE site is indicated correctly... Else die and print error
	$site !~ /\^/ and die "\nPROGRAM KILLED: Input recognition site $site missing \"^\" indicating cleavage site! \n\nFor example, the recognition site for PstI is:\n5'-CTGCA|G-3'\n3'-G|ACGTC-5'\n\nThis should be indicated 5' to 3' as:\n\"-r CTGCA^G...\"\n\nPlease call $0 -h for more help with arguments.\n\n\n";  
	#Capture fill RE site with degeneracy for use later..
	$temp = $site; 
	$temp =~ s/\^//g; 
	push (@re_plain, $temp); 
	$re_num++;
	#Check if site contains degenerate bases
	if ($site =~ /[^ACGT\^]/i){
		my @temp_re = $site; 
		print "Degenerate site $site";
		while (grep{ /[^ACGT\^]/i } @temp_re){
			for (my $t=0; $t < scalar(@temp_re); $t++){ 				
				foreach my $key (%iupac){
					if ($temp_re[$t] =~ /$key/i){ 
						my @temp_values;
						#print "Site contains $key!\n"; 
						foreach my $value (@{$iupac{$key}}){ 
							$re_pattern = $temp_re[$t];
							$re_pattern =~ s/$key/$value/i; 
							push (@temp_values, $re_pattern); 
						}
						splice (@temp_re, $t, 1, @temp_values);
					}
				}
			}
		}
		print " resolves to @temp_re\n\n"; 
		foreach my $new_re (@temp_re){ 
			$new_re =~ /([A-Z]*)\^([A-Z]*)/i;
			$re_full = $1 . $2; 
			push (@re_full, $re_full); 
			$re_pattern = $new_re; 
			$re_pattern =~ s/\^/$re_num\^$re_num/gi;
			push (@re_pattern, $re_pattern);
		}
		
		#$re_num++; 
	}else{
		#Add to array of site paterns excluding cut site marker
		#print "Input: $_\n";
		$site =~ /([A-Z]*)\^([A-Z]*)/i;
		$re_full = $1 . $2; 
		#print "Re_full: $re_full\n"; 
		push (@re_full, $re_full); 
		
		#Add to array of patterns to insert in place of above
		$re_pattern = $site; 
		$re_pattern =~ s/\^/$re_num\^$re_num/gi;
		push (@re_pattern, $re_pattern);
	
		#Increment flanking number 
		
	} 
}

print "...\n\n";
#print "RE patterns are: @re_pattern\n\n"; 

########################################################################

#Glob the input file(s)
#@input = glob "@input";

#Change directory...  
#chdir "$dirpath";

########################################################################

#In silico digest based on pattern replacement

print "\"Digesting\" genome at restriction sites: @re_full...\n\n";

	print "Reading $input... "; 
      
    #Open each of the files 
 	my @seqs;   
	$count = 0; 
	$seqs[$count]="";
	open ( FILE, "$input" ) || die "Warning: Can't open file $input!\n\n";
	    while (<FILE>){
			chomp;
		    if ($_ =~ /\>/g){   
		            $count++; 
			    $seqs[$count]=""; 
			    next;
		    }
	     
			#Store all sequences in array; different value for each 
			#print "$_\n"; 
			$seqs[$count] .= $_;  #Read file line by line and build sequence as a string       
	    }
	    
        print "$count contig\(s\) found\n";
	    close FILE;
	
	foreach (@seqs){ 
	#print "$_\n"; 
		#Make sure entry contains sequences...
		next unless $_ =~ /[A-Za-z]+/gi; 
		#Add flanking delimiters to restriction sites
		foreach (my $r=0; $r < scalar(@re_full); $r++){ 
		    #print "Replacing $re_full[$r] with $re_pattern[$r]\n"; 
		    my $full = $re_full[$r];
		    my $pattern = $re_pattern[$r];
		    $_ =~ s/$full/$pattern/gi;
		}
		   push @frags, split(/\^/, $_); #Capture all fragments in an array
	}



print "\n...\n\n"; 
 
########################################################################

#Loop through frags and build hash of hashes of frequencies for each length
#Separate based on flanking identifiers inserted with pattern replacement earlier

print "Calculating frequencies of recovered fragment lengths...\n\n"; 
#print "RE_NUM is $re_num\n"; 
#
my %fragseq; 
my $seq; 
foreach my $frag (@frags){
#print "$frag\n";
#print "Fragment: $frag\n"; 
	for (my $q=1; $q <= $re_num; $q++){ 
		for (my $w=1; $w <= $re_num; $w++){ 
			#print "Q is $q; W is $w\n";
			if ($frag =~ /^$q([A-Za-z]*?)$w$/){ 
				$len = length($1); 
				$seq = $1; 
				$key = $q . $w; 
				if ($total{$key}){ 
					$key = $q . $w; 
				}else{ 
					$key = $w . $q;
				}
				$total{$key}{$len}++;

				#If fasta files to be printed.... 
				if ($fasta==1){ 
				    push @{$fragseq{$key}}, $seq; 
				} 
			}
		}
	}  
	if ($frag !~ /^\d.*?\d$/){ 
	 
		$len = length($frag); 
		$count = $frag =~ tr/[0-9]+//; #Capture length excluding cut site IDs
			#print "Frag with $count missing flanking sites: $frag\n"; 
		$len = $len - $count; 
		$key="0";
		$total{$key}{$len}++;
		
		if ($fasta==1){ 
		    $frag =~ /^\d*([A-Za-z]*?)\d*$/; 
		    push @{$fragseq{$key}}, $1; 
		}
	}
}

#foreach my $key (keys %total){ 
#	print "$key\n"; 
#}
#print "\n"; 

########################################################################

#Check if user wants fragments printed to FASTA file

#Print all to one file, or in separate files for each fragment type?? 
#Also maybe add a size-selection mechanism for conserving file space

if ($fasta == 1){ 
my $header; 
my @names; 
my $templen; 
    print "Writing recovered fragments to FASTA file(s)...\n\n";  
    #Foreach fragment type: 
    foreach my $key (keys %total){
   	$header = "";
   	@names = split(//, $key); #Split key to get enzyme IDs
    
	#Build header and filenames 
    	if ($names[0] == 0){
    	    $header .= "Missing_sites";
    	}else{
            $header .= $re_plain[$names[0]-1] . "-" . $re_plain[$names[1]-1];
    	}
        #Open file for writing  
        open (FAS, ">$outname.$header.fasta") || die "Warning: Can't open $outname.$header.fasta!\n\n";
        
 	foreach my $seq (@{$fragseq{$key}}){
	    $templen = length($seq);
	    print FAS ">"."$header"."_Length_$templen\n";
	    print FAS "$seq"."\n";
	} 
    close FAS; 
    }
}

#Release @frags from memory
undef(@frags);
undef(%fragseq); 

########################################################################


#Create output fragment table(s)
my @temp1; 
my @temp2; 
my @merged; 
my %seen; 
$count = 0; 

my $tables = 0;

#If user chose to write to a single table... 
if ($tables == 0){ 
	
	print "Writing merged fragment table...\n\n"; 
		
	#Build list of fragment lengths to include
	#Pairwise merging of the key lists in 2D %total	
	foreach my $key (keys %total){ 

		if ($count == 0){ 
		 		
			#make array of keys from inner hash $total{$key}
			@temp1 = keys %{$total{$key}}; 
			$count++; 
			
		}else{			
			#make array of keys from *next* inner hash $total{$key}
			@temp2 = keys %{$total{$key}};
			@seen{@temp1} = ();
			
			#Make new merged array of @temp1 and nonredundant values from @temp2
			@merged = (@temp1, grep{!exists $seen{$_}} @temp2); 
			#Set @temp1 to @merged array; next iteration merge with next @temp2
			
			@temp1 = @merged; 
		}
		 
		@merged = sort {$a <=> $b} @merged; #Sort @merged numerically 
		
	}	
	
	open (OUT, ">$outname.tsv") || die "Warning: Can't open $outname.tsv: $!\n";
		
		my $header = "Fragment_length\t"; 
		my @names; 
		my %contents; 
		my $value; 
		my $sum; 
		foreach my $key (sort {$a <=> $b} keys %total){
			@names = split(//, $key); #Split key to get enzyme IDs
			#Build header based on enzyme IDs
			if ($names[0] == 0){ 
				$header .= "Missing_sites\t"; 
				
			}else{ 
				$header .= $re_plain[$names[0]-1] . "-" . $re_plain[$names[1]-1] . "\t";
			}
			foreach my $len (@merged){ 
				if ($total{$key}{$len}){ 
					$value = $total{$key}{$len};
					$contents{$len} or $contents{$len} = ""; 
					$contents{$len} .=	"$value\t";  
				}else{ 
					$value = 0;
					$contents{$len} or $contents{$len} = ""; 
					$contents{$len} .=	"$value\t"; 
				}
			}		
		}
		$header .= "Sum"; 
		print OUT "$header\n"; 
		foreach my $key (sort {$a <=> $b} keys %contents){ 
			$sum = 0; 
			$sum += $_ for split(/\s/, $contents{$key}); 			
			print OUT "$key\t" . "$contents{$key}" . "$sum\n"; 
		}
		
		close OUT;
		print "Created $outname.tsv\n\n"; 
	
}else{ 
	
	#If user chose to write individual tsv files... 
	print "Writing individual fragment tables...\n\n"; 
	my $header;
	my @names; 
	my $templen; 
	#Foreach fragment type: 
	foreach my $key (keys %total){ 
		$header = "";
		@names = split(//, $key); #Split key to get enzyme IDs
		#Build header and filenames 
		if ($names[0] == 0){ 
			$header .= "Missing_sites"; 
		}else{ 
			$header .= $re_plain[$names[0]-1] . "-" . $re_plain[$names[1]-1];
		}
		
		#Open file for writing	
 
		open ( OUT, ">$outname.$header.tsv" ) || die "Warning: Can't open $outname.$header.tsv!\n\n";

		#Print header to tsv file
		print OUT "Fragment_length\t$header\n"; 
		
		#Print fragment sizes and frequencies to table... 
		foreach my $size (sort {$a <=> $b} keys %{$total{$key}}){ 
			print OUT "\n$size\t$total{$key}{$size}";
		}
		close OUT;
		print "Created $outname.$header.tsv\n"; 
	}
	print "\n"; 
}


print "...\n\nIn silico digest complete\n\n"; 
exit;  

#################################################################################
##################################SUBROUTINES####################################
#################################################################################

#Subroutine 



#Subroutine to parse command-line arguments
sub parseArgs{

my($args) = @_;
my %opts = %$args;

my $usage="\n
Script: fragmatic.pl

Primary Author: Tyler K. Chafin - tkchafin\@uark.edu

Last Modified: 29 Jul 16

This script was created to simulate the digestion/ fragmentation of genomic DNA sequences with restriction enzymes as would be used in library preparation for a RAD-type reduced-representation genomic approach. It currently supposrts an indeterminate number of restriction enzymes, as well as degeneracy in restriction sites. 


USAGE: $0 -i /path/to/fasta_file( -r \"G^AATTC C^CGG...\" [-o /home/out][-t][-f][-m][-p]  


Mandatory Arguments
	-i	- String. Path to input FASTA file containing your genome [E.g. /home/genome.fasta] 
	-r	- String. List restriction sites in quotation marks [Usage: -r \"G^AATTC C^CGG ...\"] 

Optional Arguments
	-o	- String. Path and prefix for output files [Usage: /path/to/radsim]
	-f	- Bool. Toggle on to print fragments to FASTA files
	-h	- Bool. Display this help message 
	
fragmatic.pl will create a .tsv file containing fragment lengths for recovered loci, with the number of occurances of each length.  

";

my $i = $opts{i} or die "Input not specified!\n\n$usage\n";
my $r = $opts{r} or die "RE sites not specified!\n\n$usage\n";
my $o = $opts{o} || "sim";
my $f = $opts{f} || 0;
my $h = $opts{h} || 0;

$h == 1 and die "$usage\n";

return($i, $r, $o, $f);
	
}
