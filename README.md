# fragmatic
Simple program for in silico restriction digest of genomic sequences, to simulate RAD-family NGS library prep methods

Author: Tyler K. Chafin
Contact: tkchafin@uark.edu


##-------------------------

###Installation

This software requires Perl5 with standard modules to be installed. It assumes the perl binary to be present at /usr/bin/perl. If this is not the location of the perl binary, please edit the first line (following "#!") in fragmatic.pl accordingly. 


##-------------------------

###Usage
The help menu can be displayed using the "-h" argument:

     ./fragmatic.pl -h

The mandatory arguments are -i [input FASTA file] and -r [restriction sites]. Restriction sites are provided as a list surrounded by quotation marks, and can be any number of sites and include degenerate bases. Restriction sites containing degenerate bases will be expanded (e.g. the site for BsrFI (RCCGGY) will be expanded to ACCGGC, ACCGGT, GCCGGC, GCCGGT). You must include a caret symbol (^) indicating the cleavage location in the restriction site. To call the program using EcoRI and MspI, on an input file called "genome.fasta" located in the users home directory: 

    ./fragmatic.pl -i $HOME/genome.fasta -r "G^AATTC C^CGG"

Additional options are: 

-o: Specify a path and prefic for output files. The default is to use prefix "sim" in the current working directory. [e.g. -o /scratch/user/ddradSim]
-f: A boolean value which, when toggled "on", prints all fragments to FASTA files for further analysis. [Usage: -f]

##-------------------------

###Input file
The input file should be in standard FASTA format:

    >Scaffold1
    AGTAGATGTCGATCGATCGTACGTCG...

Contigs can be separated with their own headers, and case is not sensitive. If you have a genome assembly broken into multiple files (for each chromosome, etc), simply concatenate them together using the bash cat command:

    cat *.fasta > all.fasta

Please also note that a highly fragmentary assembly will result in many spurious fragment boundaries, or fragments which are terminated by the end of a contig rather than a restriction site. These will be tracked by the program as "Missing_sites" and can potentially obscure the fragment distribution. You can look at the frequency of these sites in the output files. One strategy to estimate the variance in locus predictions due to a highly fragmentary contig-level assembly would be to randomly concatenate the contigs across many replicates, and to perform the simulations on these replicate pseudo-assemblies. You could then estimate variance in the number of loci inferred for a given size range. I will include some scripts in the future which will help with this.


##-------------------------

###Outputs 
The standard output file is a tab-delimited table of integers representing the number of fragments of each type for a given fragment size. Using the example above, this table might look like this:

Fragment_length   | Missing_sites  |   GAATTC-GAATTC |   GAATTC-CCGG |CCGG-CCGG  |  Sum
------------------|----------------|-----------------|---------------|-----------|--------
10                |    0           |     0           |      23       |   373     |  396
11                |    0           |     1           |      63       |   384     |  448
...               |                |                 |               |           |
2846              |    0           |     0           |       0       |     1     |   1

You can also output all fragments to fasta files, which will be named according to the end-termination of each fragment (e.g. $prefix.CCGG-CCGG.fasta will contain fragments which were flanked on both ends by an MspI restriction site).


##-------------------------

###Example
A randomly generated example sequence is located in the example/ directory, as well as some outputs. To regenerate the outputs, run the fragmatic script like so, assuming that you have installed fragmatic in your home directory:

    cd $HOME/fragmatic
    ./fragmatic.pl -i example/test.fasta -o example/sim -r "G^AATTC C^CGG" -f

There is also an R script included in the scripts/ directory which can generate some plots:

    Rscript scripts/frag_plots.R example/sim.tsv

This will produce a PDF file containing several impulse plots of the resulting fragment distributions. Use these plots to qualitatively scan potential size selection ranges, trying to avoid any range likely containing a highly locus within a highly repetitive region. You can also calculate the number of loci within a given size range, and the relative proportion of sequenceable fragments, by just loading the .tsv file into excel and using the SUM function on your range and column(s) of interest. 


##-------------------------

###License
This software is provided for free: you can redistribute it and/or modify it under the terms of the GNU Public License as published by the Free Software Foundation. You should have received a copy of the GNU Public License with the software. If not, see: www.gnu.org/licenses

The author claims no liability or resposibility for the functionality of this software.
