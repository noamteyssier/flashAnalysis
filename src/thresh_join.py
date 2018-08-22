#!/usr/bin/env python3

import argparse
import sys
import difflib
import pickle

def fasta_gen(ifn):
    """generates each seq of fasta"""
    with open(ifn, 'r') as f:
        while True:
            try:
                yield [next(f)[:-1] for i in range(2)] # remove \n from information
            except StopIteration:
                break
def pair_gen(ifn):
    """generates each pair of fasta"""
    with open(ifn, 'r') as f:
        while True:
            try:
                yield [next(f)[:-1] for i in range(4)] # remove \n from information
            except StopIteration:
                break
def get_uniq(ifn):
    """
    organize a dictionary of {window : direction : sequence : count}
    organize a dictionary of {window : count}
    """
    uniq = dict() # { window : direction : sequence : count}
    uniq_count = dict() # {window : count}

    for record in fasta_gen(ifn):
        uid_full, seq = record
        uid, direction, window = uid_full.split('_')

        # add to unique window
        if window not in uniq:
            uniq[window] = dict()
            uniq_count[window] = 0

        # add to direction
        if direction not in uniq[window]:
            uniq[window][direction] = dict()

        # add sequence
        if seq not in uniq[window][direction]:
            uniq[window][direction][seq] = 0

        # increment values
        uniq[window][direction][seq] += 1
        uniq_count[window] += 1


    return (uniq, uniq_count)
def passThreshold(seq_count, winTotal, threshold):
    """ if the frequency of the sequence in the total reads passes the threshold return True"""
    frequency = float(seq_count) / winTotal * 100
    return True if frequency >= threshold else False
def hq_pair_gen(ifn, uniq, uniq_count, threshold):
    """
    calculates frequency of occurrence for each sequence in window and only
    yield high quality pairs where frequency passes threshold
    """

    for full_uid_forward, seq_forward ,full_uid_reverse, seq_reverse in pair_gen(ifn):

        # split on underscores for record information
        uid_forward, direction_forward, window_forward = full_uid_forward.split('_')
        uid_reverse, direction_reverse, window_reverse = full_uid_reverse.split('_')

        # increment counts of unique counts within window and direction
        n_fwd = uniq[window_forward][direction_forward][seq_forward]
        n_rev = uniq[window_reverse][direction_reverse][seq_reverse]

        # pull the total of high fidelity unique reads from window
        winTotal = uniq_count[window_forward]

        # if both reads pass threshold then yield
        if passThreshold(n_fwd, winTotal, threshold) and passThreshold(n_rev, winTotal, threshold):
            yield [window_forward, uid_forward, seq_forward, seq_reverse]

    return




def join_pairs(hq_pairs):
    """
    read in high quality pairs and join
    """

    memory = dict() # {seq_forward : seq_reverse : processed_overlap}
    uniq = dict()   # {window : sequence : count}


    for window, uid, seq_forward, seq_reverse in hq_pairs:

        # add new seq_forward to memory
        if seq_forward not in memory:
            memory[seq_forward] = dict()

        # add new seq_reverse pair to memory and find largest substring
        if seq_reverse not in memory[seq_forward]:
            memory[seq_forward][seq_reverse] = get_overlap(seq_forward, seq_reverse)

        # pull substring from forward and reverse pair
        substring = memory[seq_forward][seq_reverse]
        print(substring)


        try:
            seq_forward, seq_reverse, substring = find_direction(seq_forward, seq_reverse, substring)
            seq = join_seqs(seq_forward, seq_reverse, substring)

            if window not in uniq:
                uniq[window] = dict()

            if seq not in uniq[window]:
                uniq[window][seq] = 0
            uniq[window][seq] += 1

        except TypeError:
            print('misjoin %s' %window)
            # index is equal (no fix in mind yet)
            pass
    return uniq

def get_overlap(seq_forward, seq_reverse):
  """return the largest region of overlap between two seqs"""
  s = difflib.SequenceMatcher(None, seq_forward, seq_reverse)
  pos_a, pos_b, size = s.find_longest_match(0, len(seq_forward), 0, len(seq_reverse))
  return seq_forward[pos_a:pos_a+size]

def find_direction(seq_forward, seq_reverse, substring):
    """return strings in proper orientation"""
    if seq_forward.index(substring) < seq_reverse.index(substring):
        return [seq_reverse, seq_forward, substring]
    elif seq_forward.index(substring) > seq_reverse.index(substring):
        return [seq_forward, seq_reverse, substring]
    else:
        return False
def join_seqs(seq_forward, seq_reverse, overlap):
    """seq_forward until overlap, overlap, seq_reverse past overlap"""
    index1 = seq_forward.index(overlap) # position of overlap in seq_forward
    index2 = seq_reverse.index(overlap) # position of overlap in seq_reverse
    return (seq_forward[:index1] + seq_reverse[index2:])
def print_overlap(seq_forward, seq_reverse, substring):
    """format print overlap"""
    print(seq_forward)
    print('.' * seq_forward.index(substring) + substring)
    print('.' * seq_forward.index(substring) + seq_reverse)
def print_uniq(uniq):
    """print : count, sequence, direction, window"""
    for w in uniq:
        for d in uniq[w]:
            for s, c in uniq[w][d].items():
                print(c, s, d, w)
            # break
        # break
def pickle_dictionary(uniq, ofn):
    """dump uniq dictionary to file"""
    with open(ofn+'.pkl', 'wb') as f:
        pickle.dump(uniq, f, pickle.HIGHEST_PROTOCOL)

def main():
    parser = argparse.ArgumentParser(description='determine the expected snp variant for illumina seqs that are referenced to a diversity window')
    parser.add_argument('-i', '--input', help='filename of high fidelity reads to join and sort', required=True)
    parser.add_argument('-t', '--threshold', help='exclude pairs with a frequency lower than this number (default = 0.1 percent)')
    parser.add_argument('-o', '--output', help='filename to print pickled unique dictionary to', required=True)
    args = parser.parse_args()

    threshold = 0.1
    if args.threshold:
        threshold = float(args.threshold)

    uniq, uniq_count = get_uniq(args.input)
    hq_pairs = hq_pair_gen(args.input, uniq, uniq_count, threshold)
    hq_uniq_window = join_pairs(hq_pairs)

    # pickle dictionary for upstream processing
    # pickle_dictionary(hq_uniq_window, args.output)





if __name__ == '__main__':
    main()
