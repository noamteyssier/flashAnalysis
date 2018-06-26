#!/bin/bash

# assign variables
directory=$1
sample_name=$2
gen=$directory"generated_files/"

# if there is already a generated_files remove it
if [ -d $gen ] ; then
  rm -r $gen
fi

# create the generated files directory within the sample directory
mkdir $gen

# convert fastq to fasta
python fastq_to_fasta.py -d $directory -o $gen$sample_name

# grep every read with a target sequence from the fasta
./grep_reads.bash "target_seqs.fasta" $gen$sample_name".fasta" $gen"grepped_reads.fasta"

# gather high fidelity statistics and create high fidelity fasta
python hf_stats.py -g $gen"grepped_reads.fasta" -t "target_seqs.fasta" -o $gen"high_fidelity" -c $gen$sample_name".fasta"

# join high fidelity reads
python thresh_join.py -i $gen"high_fidelity.fasta" -t 0.1 -o $gen$sample_name

# gather haplotype statistics
python haplotype_pickle.py -p $gen$sample_name".pkl" -d "diversity_windows.fasta" -s "targetSnps.bed" -o $gen$sample_name


# done
echo '***Flash Analysis Complete on Sample:' $sample_name
