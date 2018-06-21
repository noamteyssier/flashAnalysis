import argparse, argcomplete, sys, os, multiprocess, time, difflib, Queue, pickle
from Bio import SeqIO
from multiprocess import Process, Queue, Pool, Manager


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
    """organize a dictionary of window -> direction -> sequence -> count"""
    uniq = dict() # { window : direction : sequence : count}
    uniq_count = dict() # {window : count}
    for pair in fasta_gen(ifn):
        uid_full, seq = [p for p in pair]
        uid, direction, window = uid_full.split('_')
        if window not in uniq:
            uniq[window] = dict()
            uniq_count[window] = 0
        if direction not in uniq[window]:
            uniq[window][direction] = dict()
        if seq not in uniq[window][direction]:
            uniq[window][direction][seq] = 0
        uniq[window][direction][seq] += 1
        uniq_count[window] += 1
    return (uniq, uniq_count)
def get_uniq_wbw(ifn):
    """uniq dictionary generator per window"""
    current_window = None
    uniq = dict()
    for u, seq in fasta_gen(ifn):
        uid, direction, window = u.split('_')
        if current_window == None: # first pass
            current_window = window
        elif current_window != window: # last pass
            # assign values
            to_return = (current_window, uniq)
            # reset values
            uniq = dict()
            current_window = window
            yield to_return
        if current_window == window: # main script
            if direction not in uniq:
                uniq[direction] = dict()
            if seq not in uniq[direction]:
                uniq[direction][seq] = 0
            uniq[direction][seq] += 1
def percentage(seq_count, window_total):
    """return the frequency of the sequence in the window specific reads as a percentage"""
    return (float(seq_count) / window_total * 100)
def hq_pair_gen(ifn, uniq, u_count, threshold):
    """generator of high quality pairs"""
    """hq_pairs = both sequences are above the treshold cutoff of percentage of total reads in window"""
    for u1,s1,u2,s2 in pair_gen(ifn):
        i1,d1,w1 = u1.split('_')
        i2,d2,w2 = u2.split('_')
        count_1 = uniq[w1][d1][s1]
        count_2 = uniq[w2][d2][s2]
        win_total = u_count[w1]
        if (percentage(count_1, win_total) >= threshold) and (percentage(count_2, win_total) >= threshold):
            yield u1, s1, u2, s2
            # break
def join_pairs(hq_pairs):
    memory = dict()
    uniq = dict()
    for u1,s1,u2,s2 in hq_pairs:
        pid, direction, window = u1.split('_')
        # uid = '_'.join([pid, window])
        if s1 not in memory:
            memory[s1] = dict()
        if s2 not in memory[s1]:
            memory[s1][s2] = get_overlap(s1, s2)
        ss = memory[s1][s2]
        try:
            s1, s2, ss = find_direction(s1, s2, ss)
            seq = join_seqs(s1, s2, ss)
            if window not in uniq:
                uniq[window] = dict()
            if seq not in uniq[window]:
                uniq[window][seq] = 0
            uniq[window][seq] += 1
        except TypeError:
            # index is equal (no fix in mind yet)
            pass
    return uniq


def get_overlap(s1, s2):
  """return the largest region of overlap between two seqs"""
  s = difflib.SequenceMatcher(None, s1, s2)
  pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
  return s1[pos_a:pos_a+size]
def find_direction(s1, s2, ss):
    """return strings in proper orientation"""
    if s1.index(ss) < s2.index(ss):
        return [s2, s1, ss]
    elif s1.index(ss) > s2.index(ss):
        return [s1, s2, ss]
    else:
        return False
def join_seqs(s1, s2, overlap):
    """s1 until overlap, overlap, s2 past overlap"""
    index1 = s1.index(overlap) # position of overlap in s1
    index2 = s2.index(overlap) # position of overlap in s2
    return (s1[:index1] + s2[index2:])
def print_overlap(s1, s2, ss):
    """format print overlap"""
    print s1
    print '.' * s1.index(ss) + ss
    print '.' * s1.index(ss) + s2
def print_uniq(uniq):
    """print : count, sequence, direction, window"""
    for w in uniq:
        for d in uniq[w]:
            for s, c in uniq[w][d].items():
                print c, s, d, w
            # break
        # break
def pickle_dictionary(uniq, ofn):
    """dump uniq dictionary to file"""
    with open(ofn+'.pkl', 'wb') as f:
        pickle.dump(uniq, f, pickle.HIGHEST_PROTOCOL)

def main():
    parser = argparse.ArgumentParser(description='determine the expected snp variant for illumina seqs that are referenced to a diversity window')
    parser.add_argument('-i', '--input', help='filename of high fidelity reads to join and sort', required=True)
    parser.add_argument('-t', '--threshold', help='exclude pairs with a percentage lower than this number (default = 0.1 percent)')
    parser.add_argument('-o', '--output', help='filename to print pickled unique dictionary to', required=True)
    args = parser.parse_args()

    threshold = 0.1
    if args.threshold:
        threshold = float(args.threshold)

    uniq, uniq_count = get_uniq(args.input)
    hq_pairs = hq_pair_gen(args.input, uniq, uniq_count, threshold)
    hq_uniq_window = join_pairs(hq_pairs)

    # pickle dictionary for upstream processing
    pickle_dictionary(hq_uniq_window, args.output)





if __name__ == '__main__':
    main()
