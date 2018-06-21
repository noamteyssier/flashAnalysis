import os, sys, argparse, argcomplete, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

# NEED TO FIX:
# PRINT AS REVERSE COMPLEMENT IF TARGET SEQUENCE IS REVERSE

class seqio:
    def __init__(self, id, seq):
        self.id = id[1:-1]
        self.seq = seq[:-1]
        self.pair_bool = False
        self.target = None
    def add_target(self, target):
        self.target = target
    def switch_bool(self):
        self.pair_bool = True
    def write_fasta(self, window, fn):
        fn.write('>' + self.id + '_' + window.id + '\n')
        if '_forward' in self.target.id:
            fn.write(self.seq + '\n')
        else:
            quick_seq = Seq(self.seq)
            fn.write(str(quick_seq.reverse_complement()) + '\n')
class Diversity_Window:
    windows = list()
    dictionary = dict()
    def __init__(self, id):
        self.id = id
        self.targets = list()
        self.hf_count = 0
        self.left = None
        self.right = None
        Diversity_Window.windows.append(self)
    def add_target(self, target):
        """add to list of targets"""
        self.targets.append(target)
        if 'start_forward' in target.id:
            self.left = target
        elif 'end_reverse' in target.id:
            self.right = target
    def increment_hf(self):
        self.hf_count += 1
    def left_cuts(self):
        return self.left.total_reads
    def right_cuts(self):
        return self.right.total_reads
    def total_cuts(self):
        return self.left_cuts() + self.right_cuts()
    def print_stats(self):
        stats = [self.total_cuts(), self.left_cuts(), self.right_cuts(), self.hf_count, self.id]
        return [str(s) for s in stats]
class Target:
    targets = list()
    lookup = dict()
    lookup_full_seq = dict()
    def __init__(self, seqio):
        """add self to target list and into dictionary for lookup"""
        self.seqio = seqio
        self.id = seqio.id
        self.seq = str(seqio.seq)
        self.window = None

        # master list of reads
        self.pairs = list()

        # total counts
        self.total_reads = 0
        self.z_count = 0
        self.o_count = 0
        self.high_fidelity = 0

        Target.targets.append(self)
        Target.lookup[self.seq] = self
    def assign_window(self, window):
        """assign window to target for easy link"""
        self.window = window
    def add_read(self, seq, pos=0):
        """add read to list of reads"""
        self.total_reads += 1
        if pos == 0:
            self.z_count += 1
        if pos == 1:
            self.o_count += 1
        Target.lookup_full_seq[seq] = self
    def add_pair(self, pair):
        self.pairs.append(pair)
        self.high_fidelity += 1

def parse_targets(tfn):
    """create target objects and diversity window objects"""
    """assign targets to their windows"""
    for r in SeqIO.parse(tfn, 'fasta'):
        if ('start_reverse' not in r.id) and ('end_forward' not in r.id) and ('alternate' not in r.id):
            t = Target(r)
            w = r.id.split('_')[0]
            if w not in Diversity_Window.dictionary:
                Diversity_Window.dictionary[w] = Diversity_Window(w)
            Diversity_Window.dictionary[w].add_target(t)
            t.assign_window(Diversity_Window.dictionary[w])
def fasta_reader(fn):
    """generate two lines of fasta"""
    with open(fn, 'r') as f:
        while True:
            try:
                id, seq = [next(f)[:-1] for i in range(2)]
                yield id,seq
            except StopIteration:
                break
def assign_grepped(gfn):
    """assign generated reads to a target"""
    pairs = dict() # {window : uid : seq}
    for id, seq in fasta_reader(gfn):
        z_pos = seq[:20]
        o_pos = seq[1:21]
        try: # target sequence found in position 0
            t = Target.lookup[z_pos]
            position = 0
        except KeyError:
            try: # target sequence found in position 1
                t = Target.lookup[o_pos]
                position = 1
            except KeyError: # no target sequence is found
                continue

        # split uid
        uid, direction = id.split('_')

        # add read to target sequence to gather stats
        t.add_read(seq, pos=position)

        # organize pairs in nested dictionary
        if t.window.id not in pairs:
            pairs[t.window.id] = dict()
        if uid not in pairs[t.window.id]:
            pairs[t.window.id][uid] = list()

        # append direction to sequence for later retrieval
        pairs[t.window.id][uid].append(direction + '_' + seq)

        # increment high fidelity read count by two when pair is found
        if len(pairs[t.window.id][uid]) == 2:
            t.window.increment_hf()

    return pairs
def total_reads(cat):
    """grep number of sequences and fasta and return int"""
    return int(subprocess.check_output(['fgrep ">" -c ' + cat], shell=True))
def print_stats(output, total):
    with open(output+'.csv', 'w+') as f:
        f.write('total_reads, total_cuts, left_cuts, right_cuts, high_fidelity_pairs, name\n')
        for w in Diversity_Window.windows:
            f.write(str(total) + ',' + ','.join(w.print_stats()) + '\n')
def hf_fasta(output, pairs):
    with open(output+'.fasta', 'w+') as f:
        for window in pairs:
            for uid in pairs[window]:
                if len(pairs[window][uid]) == 2: # only print high fidelity reads
                    for s in pairs[window][uid]:
                        direction, seq = s.split('_') # separate direction from seq (added in assigned_grepped)
                        header = [uid, direction, window]
                        if '_reverse' in Target.lookup_full_seq[seq].id: # write reverse complement
                            quick_seq = Seq(seq)
                            seq = str(quick_seq.reverse_complement())
                        f.write('_'.join(header) + '\n' + seq + '\n')


def main():
    parser = argparse.ArgumentParser(description='create high fidelity statistics')
    parser.add_argument('-g', '--grepped', help='filename of grepped reads', required=True)
    parser.add_argument('-t', '--targets', help='filename of target sequences', required=True)
    parser.add_argument('-o', '--output', help='filename of output files (will make csv and fasta)', required=True)
    parser.add_argument('-c', '--cat', help='filename of concatenated reads to generate total read count', required=True)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # assign args
    grepped = args.grepped
    targets = args.targets
    output = args.output
    cat = args.cat

    # methods
    parse_targets(targets)
    pairs = assign_grepped(grepped)
    total = total_reads(cat)

    # printing methods
    print_stats(output, total)
    hf_fasta(output, pairs)



if __name__ == '__main__':
    main()
