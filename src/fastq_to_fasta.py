#!/usr/bin/env python

import Bio
import sys
import os
import argparse
import re
import gzip
from Bio import SeqIO


def get_filenames(directory):
    """return the filenames of the fast(q/a)s or complain directory doesnt exist"""

    try:
        # match contents of directory to regex and assemble list
        fastq_list = [directory + i for i in os.listdir(directory) if re.match('[a-zA-Z0-9-_.]+.(fastq|fq)', i)]

        # only return if two fastqs are in directory
        if len(fastq_list) == 2:
            return fastq_list
        else:
            sys.exit('Must only be 2 fastq/fq samples in directory')

    except OSError:
        sys.exit('Cannot Find Given Directory <%s>' %directory)

def q2a(filename, direction, f):
    """convert fastq to fasta"""

    # check for gzip extension and open accordingly
    fastq = gzip.open(filename, 'r') if '.gz' in filename \
        else open(filename, 'r')

    while True:
        try:
            uid, seq, plus, phred = [next(fastq) for i in range(4)]
            uid = uid.split(' ')[0].replace('@', '>')
            if direction == 'F':
                f.write(uid + '_1\n')
            elif direction == 'R':
                f.write(uid + '_2\n')
            f.write(seq)
        except StopIteration:
            break

def write_fastas(forward, reverse, output):
    f = open(output, 'w+')
    q2a(forward, 'F', f)
    q2a(reverse, 'R', f)

def main():
    parser = argparse.ArgumentParser(description='convert two illumina fastq files into three fastas (forward, reverse, and both)')
    parser.add_argument('-d', '--directory', help='directory containing fastq files (can only have two)', required=True)
    parser.add_argument('-o', '--output', help='output filename for converted file', required=True)
    args = parser.parse_args()



    # perform fastq -> fasta
    fn1, fn2 = get_filenames(args.directory)
    write_fastas(fn1, fn2, args.output)

if __name__ == '__main__':
    main()
