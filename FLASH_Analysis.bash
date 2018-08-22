#!/usr/bin/env bash

# assign variables
directory=$1
sample_name=$2
gen=$directory"generated_files/"

# if there is already a generated_files remove it
if [ -d $gen ] ; then
  rm -r $gen
fi
#
# create the generated files directory within the sample directory
mkdir $gen

# convert fastq to fasta
python src/fastq_to_fasta.py -d $directory -o $gen$sample_name".fa"

# # gather high fidelity statistics and create high fidelity fasta
python3 src/hf_stats.py \
  -t "target_seqs.tab.txt"  \
  -i $gen$sample_name".fa"  \
  -b $sample_name           \
  -d $gen

# # join high fidelity reads
# python thresh_join.py -i $gen"high_fidelity.fasta" -t 0.1 -o $gen$sample_name
#
# # gather haplotype statistics
# python haplotype_pickle.py -p $gen$sample_name".pkl" -d "diversity_windows.fasta" -s "targetSnps.bed" -o $gen$sample_name
#
#
# # done
# echo '***Flash Analysis Complete on Sample:' $sample_name
