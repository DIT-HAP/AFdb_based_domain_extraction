#!/bin/bash

# This script takes the output of the PAE analysis and converts it to a domain

# Input: 1) PAE output file
# Output: 1) Domain file

inputFolder=$1
outputFolder=$2

# Create output folder
mkdir -p $outputFolder


# Get the PAE output file
PAE_files=`find $inputFolder -name "*.json"`

# Loop through each PAE output files
for PAE_file in $PAE_files; do

    # Get the domain name
    domain_name=`basename $PAE_file | cut -d "." -f 1`
    domain_name=`basename $PAE_file | sed s/.json//g`

    # Get the domain file with multiple processes
    domain_file=$outputFolder/$domain_name".domain"
    nohup python3 pae_to_domains_fix_random_seed.py --resolution 0.5 --output_file $domain_file $PAE_file &

done