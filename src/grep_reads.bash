#!/usr/bin/env bash

targetSeqs=$1   # target_seqs filename
genFiles=$2     # generated_files directory
input_fa=$3     # input fasta filename

# pull targetSeqs into tab delim and apply regex pattern
# will match target sequence if it begins at zero or one position only
# will also remove sequences with N's
cat $targetSeqs | \
  cut -f 2,3 | \
  sed 's/\t/\n/' | \
  sed 's/^/\t[ACTG]?/;s/$/[ACTG]+$/' > $genFiles"reg_targets.txt"

# grep using regular expressions and pipe to output
cat $input_fa | \
  paste - - | \
  egrep -f $genFiles"reg_targets.txt" | \
  sed 's/\t/\n/'
