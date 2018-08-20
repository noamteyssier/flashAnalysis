
t=$1        # target_seqs filename
fa=$2       # input fasta filename
o=$3        # output grep_reads filename

fgrep -f $t -B 1 --no-group-separator $fa > $o
