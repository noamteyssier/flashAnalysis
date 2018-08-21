#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess

from modules import Diversity_Window

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
            'src/grep_reads.bash %s %s %s' %(target_filename, genDir, input_filename), \
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
def write_stats(statsOut, totalReads):
    """
    Write cut statistics for each window to a tab delim file
    """
    header = [
        'totalReads',
        'totalCuts',
        'leftCuts',
        'rightCuts',
        'highFidelityCount',
        'window'
        ]

    statsOut.write('\t'.join(header)+'\n')
    for window in Diversity_Window.windows:
        writeOut = [
            totalReads,
            window.get_totalCuts(),
            window.get_leftCuts(),
            window.get_rightCuts(),
            window.get_highFidelityCount(),
            window.name
            ]
        statsOut.write('\t'.join(writeOut)+'\n')
def write_fasta(fastaOut, window, uid, forwardSeq, reverseSeq):
    """
    write a fasta entry for each sequence giving appropriate direction
    and window information
    """
    fasta_entry = [
    '_'.join(
        [uid, 'forward', window.name]
        ),
    forwardSeq,
    '_'.join(
        [uid, 'reverse', window.name]
        ),
    reverseSeq
    ]

    fastaOut.write('\n'.join(fasta_entry)+'\n')
def main():
    parser = argparse.ArgumentParser(description='create high fidelity statistics')
    parser.add_argument('-t', '--targets', help='filename of target sequences', required=True)
    parser.add_argument('-b', '--baseName', help='baseName to preappend to tab.txt and fasta output', required=True)
    parser.add_argument('-i', '--input', help='fasta to grep reads and assign to targets')
    parser.add_argument('-d', '--directory', help='generated_files directory to write to')
    args = parser.parse_args()

    # assign args
    input = args.input
    targets = args.targets
    genDir = args.directory
    baseName = args.baseName
    outName = genDir+baseName+'_hf'

    statsOut = open(outName+'.tab.txt', 'w+')
    fastaOut = open(outName+'.fa', 'w+')

    # # methods
    parse_targets(targets)
    totalReads = calculate_total(input)
    for window, uid, forwardSeq, reverseSeq in assign_grepped(input, targets, genDir):
        write_fasta(fastaOut, window, uid, forwardSeq, reverseSeq)

    write_stats(statsOut, totalReads)



if __name__ == '__main__':
    main()
