import os, sys, Bio, argparse, modules, tempfile, pickle
from Bio import SeqIO, AlignIO
from modules import * # MIP, Diversity_Window, SNP, Alignment, Fasta, Pair
from subprocess import Popen, PIPE, STDOUT


def depickle(pkl):
    """depickle an object and return"""
    return pickle.load(open(pkl, 'rb'))
def make_windows(dw_fn):
    """make window objects for diversity windows"""
    for r in SeqIO.parse(dw_fn, 'fasta'):
        d = Diversity_Window()
        d.add_seqio(r)
def make_snps(snp_fn):
    """make snp objects for each snp found in tsv"""
    with open(snp_fn, 'r') as f:
        next(f) # skip header line
        for line in f:
            s = SNP()
            s.add_line(line)
def chromosome_lookup(dw=None, snp=None):
    """return chromosome of dw given snp and vice versa"""
    if snp != None:
        return [d for d in Diversity_Window.windows if d.chrom_num == snp.chrom]
def assign_snps_to_windows():
    """assign each snp to a diversity window"""
    for s in SNP.snps:
        [d.add_snp(s) for d in chromosome_lookup(snp=s)]
def assign_uniq_to_windows(uniq):
    """will assign each uniq dictionary to respective window"""
    for u in uniq:
        Diversity_Window.lookup[u].add_uniq(uniq[u])
def join_and_assign_reads(fn):
    """parse high fidelity reads, join on overlap, and add to window"""
    for r in SeqIO.parse(fn, 'fasta'):
        uid = r.id.split('_')[0]
        if uid not in Pair.lookup:
            Pair(uid)
        if Pair.lookup[uid].add_seq(r) == True:
            window = r.id.split('_')[2]
            Diversity_Window.lookup[window].add_read(Pair.lookup[uid])
def align_reads(hf_threshold=0.01):
    """align high fidelity reads to window and add corresponding snp strings to windows"""
    for w in Diversity_Window.windows:
        if len(w.unique_reads) > 0:
            a = Alignment()
            a.hf_and_dw(w.hf_gen(threshold=hf_threshold), w)
            [w.add_snp_string(s) for s in a.get_snps_columns()]
            [w.add_haplotype(h) for h in a.get_snps_rows()]
def count_snps():
    """count snps for windows"""
    for w in Diversity_Window.windows:
        if len(w.unique_reads) > 0:
            w.count_snps()
def print_haplotypes(ofn):
    with open(ofn + '_haplo.csv', 'w+') as f:
        f.write('total_reads, percentage, haplotype_string, window\n')
        for w in Diversity_Window.windows:
            if len(w.unique_haplotypes) > 0:
                w.print_haplotypes(f)
def print_snps(ofn):
    with open(ofn + '_snps.csv', 'w+') as f:
        f.write('(chrom:pos), total reads, percent_ref, percent_nref, percent_unknown\n')
        for w in Diversity_Window.windows:
            w.print_snp_count(f)

def main():
    parser = argparse.ArgumentParser(description='determine the expected snp variant for illumina seqs that are referenced to a diversity window')
    parser.add_argument('-p', '--pickle', help='filename of pickled unique dictionary', required=True)
    parser.add_argument('-d', '--dw_fn', help='filename of Diversity Windows to align against', required=True)
    parser.add_argument('-s', '--snp_fn', help='filename of SNPS to look for in alignment', required=True)
    parser.add_argument('-o', '--output', help='filename to print snps to (will print as csv)', required=True)
    args = parser.parse_args()

    # make objects
    make_windows(args.dw_fn)
    make_snps(args.snp_fn)
    uniq = depickle(args.pickle)

    # assign objects to hierarchy
    assign_snps_to_windows()
    assign_uniq_to_windows(uniq)


    # align reads to window and count snps
    align_reads()
    count_snps()

    # print the snp ratios for all snps
    print_snps(args.output)
    print_haplotypes(args.output)


if __name__ == '__main__':
    main()
