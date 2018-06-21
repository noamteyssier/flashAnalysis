# Author : Noam Teyssier
#
# Shell Script to grep reads from fasta
# example usage:
#       sh grep_reads.bash target_seqs.fasta input/edc4.fasta grepped_reads.fasta
#./grep_reads.bash generated_files/target_seqs.fasta generated_files/edc4.fasta generated_files/grepped_reads.fasta generated_files/outside_cuts.fasta generated_files/outside_reads.fasta

t=$1        # target_seqs filename
fa=$2       # input fasta filename
o=$3        # output grep_reads filename
#outside=$4  # outside cuts filename
#cutsfn=$5   # outside grep cuts filename
fgrep -f $t -B 1 --no-group-separator $fa >$o
#echo reads written to... '          <'$o'>'
#fgrep -f $outside -B 1 --no-group-separator $fa >$cutsfn
#echo outside_cuts written to... '   <'$cutsfn'>'
