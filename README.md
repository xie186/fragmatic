# fragmatic
Simple program for in silico restriction digest of genomic sequences, to simulate RAD-family NGS library prep methods

Author: Tyler K. Chafin
Contact: tkchafin@uark.edu

###Installation

This software requires Perl5 with standard modules to be installed. It assumes the perl binary to be present at /usr/bin/perl. If this is not the location of the perl binary, please edit the first line (following "#!") in the fragmatic.pl accordingly. 

###Usage
The help menu can be displayed using the "-h" argument:
./fragmatic.pl -h

The mandatory arguments are -i [input FASTA file] and -r [restriction sites]. Restriction sites are provided as a list surrounded by quotation marks, and can be any number of sites and include degenerate bases. Restriction sites containing degenerate bases will be expanded (e.g. the site for BsrFI (RCCGGY) will be expanded to ACCGGC, ACCGGT, GCCGGC, GCCGGT). You must include a caret symbol (^) indicating the cleavage location in the restriction site. To call the program using EcoRI and MspI, on an input file called "genome.fasta" located in the users home directory: 

./fragmatic.pl -i $HOME/genome.fasta -r "G^AATTC C^CGG"

Additional options are: 

-o: Specify a path and prefic for output files. The default is to use prefix "sim" in the current working directory. [e.g. -o /scratch/user/ddradSim]
-f: A boolean value which, when toggled "on", prints all fragments to FASTA files for further analysis. [Usage: -f]

##Input file
The input file should be in standard FASTA format, like so:

>Scaffold1
AGTAGTGATGATGATGAGAGAGATA...

Contigs can be separated with their own headers, and case is not sensitive. If you have a genome assembly broken into multiple files (for each chromosome, etc), simply concatenate them together using the bash cat command:

cat *.fasta > all.fasta

Please also note that a highly fragmentary assembly will result in many spurious fragment boundaries, or fragments which are terminated by the end of a contig rather than a restriction site. These will be tracked by the program as "Missing_sites" and can potentially obscure the fragment distribution. You can look at the frequency of these sites in the output files.

##Outputs 
The standard output file is a tab-delimited table of integers representing the number of fragments of each type for a given fragment size. Using the example above, this table might look like this:

Fragment_length    Missing_sites    GAATTC-GAATTC    GAATTC-CCGG CCGG-CCGG    Sum
10                    0                0                 23          373      396
11                    0                1                 63          384      448
...
2846                  0                0                 0            1        1








###License
This software is provided for free: you can redistribute it and/or modify it under the terms of the GNU Public License as published by the Free Software Foundation. You should have received a copy of the GNU Public License with the software. If not, see: www.gnu.org/licenses

The author claims no liability or resposibility for the functionality of this software.
