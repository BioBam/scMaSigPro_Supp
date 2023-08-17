#!/bin/bash

# Define the input file
INPUT_FILE="Pubmed-summary-Trajectory-set.txt"

# Define the output file
OUTPUT_FILE="dois.txt"

# Use grep to find lines containing 'doi', then use sed to isolate the DOI.
grep -o 'doi: [^ ]*' $INPUT_FILE | sed 's/^doi: //' > $OUTPUT_FILE
