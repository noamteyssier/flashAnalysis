#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

from modules import Diversity_Window

# NEED TO FIX:
# PRINT AS REVERSE COMPLEMENT IF TARGET SEQUENCE IS REVERSE
#
# class seqio:
#     def __init__(self, id, seq):
#         self.id = id[1:-1]
#         self.seq = seq[:-1]
#         self.pair_bool = False
#         self.target = None
#     def add_target(self, target):
#         self.target = target
#     def switch_bool(self):
#         self.pair_bool = True
#     def write_fasta(self, window, fn):
#         fn.write('>' + self.id + '_' + window.id + '\n')
#         if '_forward' in self.target.id:
#             fn.write(self.seq + '\n')
#         else:
#             quick_seq = Seq(self.seq)
#             fn.write(str(quick_seq.reverse_complement()) + '\n')
# class Diversity_Window:
#     windows = list()
#     dictionary = dict()
#     def __init__(self, id):
#         self.id = id
#         self.targets = list()
#         self.hf_count = 0
#         self.left = None
#         self.right = None
#         Diversity_Window.windows.append(self)
#     def add_target(self, target):
#         """add to list of targets"""
#         self.targets.append(target)
#         if 'start_forward' in target.id:
#             self.left = target
#         elif 'end_reverse' in target.id:
#             self.right = target
#     def increment_hf(self):
#         self.hf_count += 1
#     def left_cuts(self):
#         return self.left.total_reads
#     def right_cuts(self):
#         return self.right.total_reads
#     def total_cuts(self):
#         return self.left_cuts() + self.right_cuts()
#     def print_stats(self):
#         stats = [self.total_cuts(), self.left_cuts(), self.right_cuts(), self.hf_count, self.id]
#         return [str(s) for s in stats]
# class Target:
#     targets = list()
#     lookup = dict()
#     lookup_full_seq = dict()
#     def __init__(self, seqio):
#         """add self to target list and into dictionary for lookup"""
#         self.seqio = seqio
#         self.id = seqio.id
#         self.seq = str(seqio.seq)
#         self.window = None
#
#         # master list of reads
#         self.pairs = list()
#
#         # total counts
#         self.total_reads = 0
#         self.z_count = 0
#         self.o_count = 0
#         self.high_fidelity = 0
#
#         Target.targets.append(self)
#         Target.lookup[self.seq] = self
#     def assign_window(self, window):
#         """assign window to target for easy link"""
#         self.window = window
#     def add_read(self, seq, pos=0):
#         """add read to list of reads"""
#         self.total_reads += 1
#         if pos == 0:
#             self.z_count += 1
#         if pos == 1:
#             self.o_count += 1
#         Target.lookup_full_seq[seq] = self
#     def add_pair(self, pair):
#         self.pairs.append(pair)
#         self.high_fidelity += 1

def parse_targets(target_filename):
    """
    - create Target & Diversity_Window classObjects
    - assign target sequences to appropirate windows
    """
    f = open(target_filename, 'r')
    while True:
        try:
            window, startForward, endReverse = next(f).strip('\n').split('\t')
            w = Diversity_Window(window)
            w.add_targets(startForward, endReverse)
        except StopIteration:
            return


    # for r in SeqIO.parse(tfn, 'fasta'):
    #     t = Target(r)
    #     w = r.id.split('_')[0]
    #     if w not in Diversity_Window.dictionary:
    #         Diversity_Window.dictionary[w] = Diversity_Window(w)
    #     Diversity_Window.dictionary[w].add_target(t)
    #     t.assign_window(Diversity_Window.dictionary[w])
def calculate_total(input_filename):
    wc_l = subprocess.Popen(
        'wc -l %s' %input_filename,
        shell=True, stdout=subprocess.PIPE
    )

    out,err = wc_l.communicate()
    total = out.decode('ascii').split(' ')[0]
    return total
def grepReads(input_filename, target_filename, genDir):
    """call grep with regex and pipe directly"""
    grep_out = subprocess.Popen(
            'src/new_grep_reads.bash %s %s %s' %(target_filename, genDir, input_filename), \
            stdout=subprocess.PIPE, shell=True
        )

    f = iter(grep_out.stdout.readline, 'b')
    while True:
        id, seq = [next(f).decode('ascii').strip('\n') for i in range(2)]
        if id:
            yield id,seq
        else:
            break
def windowLookup(grepSeq):
    """iterate through windows to assign sequence"""
    for window in Diversity_Window.windows:
        if window.matchToTarget(grepSeq):
            return window
    return False
def assign_grepped(input, targets, genDir):
    """
    generator that :
    - assigns generated reads to a target
    - increments values in assigned Diversity_Window
    - yields [window, uid, forwardSeq, reverseSeq]
    """
    # pairs = dict() # {window : uid : seq}

    grepSeq_to_seqID = dict() # {grepSeq : seqID}
    seqID_to_grepSeq = dict() # {seqID : grepSeq}
    seqID_to_window = dict() # {seqID : window}
    pair_uid = dict() # {uid : [seqID_1, seqID_2]}

    seqID_tally = 0

    for id, grepSeq in grepReads(input, targets, genDir):

        # give grepSeq a seqID if not seen before and make accesible for lookup
        if grepSeq not in grepSeq_to_seqID:

            # assign seqID and increment tally
            seqID = seqID_tally
            seqID_tally += 1

            # create new entries in lookup dicts
            grepSeq_to_seqID[grepSeq] = seqID
            seqID_to_grepSeq[seqID] = grepSeq

            # perform window lookup and push to seqID_to_window
            window = windowLookup(grepSeq)
            seqID_to_window[seqID] = window

        # grepSeq processed before, pull window from lookup
        else:
            seqID = grepSeq_to_seqID[grepSeq]
            window = seqID_to_window[seqID]

        if window:

            # get uid and seqID
            uid = id.split('_')[0]

            # append to pair_uid
            if uid not in pair_uid:
                pair_uid[uid] = []
            pair_uid[uid].append(seqID)

            # yield only pairs where both reads were assigned
            if len(pair_uid[uid]) == 2:
                window.increment_highFidelityCount()
                yield [window, uid] + [seqID_to_grepSeq[seqID] for seqID in pair_uid[uid]]
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
    parser.add_argument('-b', '--basename', help='basename to preappend to tab.txt and fasta output', required=True)
    parser.add_argument('-i', '--input', help='fasta to grep reads and assign to targets')
    parser.add_argument('-d', '--directory', help='generated_files directory to write to')
    args = parser.parse_args()

    # assign args
    input = args.input
    targets = args.targets
    genDir = args.directory
    basename = args.basename

    ## these arguments would get around having to ask for another argument to calculate total reads
    # run grep reads as subprocess
    # args.input = fasta to grep target sequences in
    # args.targets = tab.txt to create regex and grep with
    # args.basename = basename to preappend to output files


    # # methods
    parse_targets(targets)
    calculate_total(input)
    for assignedPairList in assign_grepped(input, targets, genDir):
        window, uid, forwardSeq, reverseSeq = assignedPairList
        break


    # # printing methods
    # print_stats(output, total)
    # hf_fasta(output, pairs)
    #


if __name__ == '__main__':
    main()
