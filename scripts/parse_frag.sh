#!/bin/bash

#Check command line arguments and assign
if [ $1 ]; then 
  file="$1"
else
  printf "\nUsage: $0 <File Name> <LocusColumn (int)> <From (int)> <To (int)>\n"
  printf "This script is means to parse the .tsv output of fragmatic.pl\n\n"
  exit 1
fi

#Check that file exists
if [ ! -e $1 ]; then 
  printf "\nFile $1 does not exist!\n"
  printf "\nUsage: $0 <File Name> <LocusColumn #>\n" 
  printf "This script is means to parse the .tsv output of fragmatic.pl\n\n"
  exit 1; 
fi;

#Prompt for column to use (if not provided)
if [ $2 ]; then
  col="$2"
else
  printf "\nNo column choice provided. Which column would you like to parse (enter #)?\n\nChoices:\n"
  head -1 $1 | awk 'BEGIN{RS="\t";num=1}{print "[" num "]: " $0; num++}'
  read -p "Column: " col
fi

#Prompt for column to use (if not provided)
if [ $3 ] && [ $4 ]; then
  from=$3
  to=$4
  if [ $to -lt $from ] || [ $to -eq $from ]; then
    printf "Error! Size selection upper limit must be greater than lower limit!\n"
    exit 1
  fi
else
  printf "\nSize selection not specified\n"
  read -p "Lower bound (int): " from
  read -p "Upper bound (int): " to
fi

name=`head -1 $file | awk -v col=$col '{print $col}'`

printf "\nInput: \n  Parsing column <$name> in file <$file>. \n  Size selection: $from - $to\n\n";
echo "-----------------------------------------------"
printf "Result: \n"
cat $file | awk '{
	if ($1 > "'"$from"'" && $1 < "'"$to"'"){ 
		loc+=$"'"$col"'"; 
		sum+=$6
	}
}END{
	print "  Number of sequence-able fragment: "loc; 
	print "  Proportion of sequence-able fragments: " loc/sum
}' 
printf "\n\n"
