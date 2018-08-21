import sys
import os
import subprocess
import tempfile
import difflib
import time
import re
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from subprocess import Popen, PIPE, STDOUT

class MIP:
    mip_names = list()
    mips = list()
    def __init__(self, id=None):
        self.id = id
        self.sequences = list()
        self.strains = list()
        self.expected_snps = list()
        self.diversity_window = None
        MIP.mips.append(self)
        if self.id != None:
            MIP.mip_names.append(self.id)
    def add_sequence_strains(self, seq, strains):
        """add sequence and strains to respective lists and preserve index"""
        self.sequences.append(seq)
        self.strains.append(strains)
    def find_sequence_given(self, strain):
        """find the index of the strain and return the sequence at that index"""
        for i in range(len(self.strains)):
            if strain in self.strains[i]:
                return self.sequences[i]
    def add_seqio(self, seqio):
        """add seqio information to self"""
        self.id = seqio.id
        self.seq = str(seqio.seq)
        self.seqio = seqio
    def add_diversity_window(self, dw):
        self.diversity_window = dw
    def parse_dw_from_id(self):
        """return the diversity window id from the mip id"""
        self.name, self.dw_name, self.direction = self.id.split(':__:')
        return self.dw_name
    def add_expected(self, dw, exp, pos=0):
        """add expected snps to class"""
        self.expected_snps = [e[pos] for e in exp]
    def csv_expected(self):
        """print csv of mipname, snp_pos, ref, nref, strain"""
        # mip name, snp position, reference variant, nonreference variant, strain variant
        dw_snps = self.diversity_window.snps
        for i in range(len(self.expected_snps)):
            info = [self.name, dw_snps[i].full_position, dw_snps[i].ref, dw_snps[i].nref, self.expected_snps[i]]
            str_info = [str(i) for i in info]
            print(','.join(str_info))
class Fasta:
    def dw_write(self, fn, dw):
        """write diversity window given a filename"""
        with open(fn, 'w+') as f:
            f.write('>' + str(dw.id) + '\n')
            f.write(str(dw.seq) + '\n')
    def mip_write(self, fn, mip, direction=False):
        """write to filename given mip and if direction matters"""
        with open(fn, 'a+') as f:
            f.write('>' + str(mip.id) + '\n')
            if direction == True:
                if 'reverse' in mip.id:
                    f.write(str(mip.seqio.seq.reverse_complement()) + '\n')
                    return
            f.write(str(mip.seq) + '\n')
class Diversity_Window:
    windows = list()
    window_ids = list()
    lookup = dict()
    def __init__(self, name):
        self.name = name
        self.id = None
        self.seq = None
        self.seqio = None
        self.chrom = None
        self.chrom_num = None
        self.left_pos = None
        self.right_pos = None

        # target information
        self.startForward = None
        self.endReverse = None

        ## counts
        self.leftCuts = 0
        self.rightCuts = 0
        self.highFidelityCount = 0
        self.total_reads = 0
        self.total_haplo = 0


        # organizers
        self.reads = dict() # {seq : count}


        self.snps = list()
        self.mips = list()
        self.unique_reads = dict() # dictionary to store unique reads and their counts
        self.unique_haplotypes = dict() # dictionary to store unique haplotypes and their counts
        self.snp_strings = list()  # collection of snp strings taken from alignment
        self.alignment_order = list() # the ids in order of alignment
        self.aln_order_weights = list() # the weights isolated from alignment ids
        Diversity_Window.windows.append(self)
    def add_seqio(self, seqio):
        self.id = seqio.id
        self.seq = str(seqio.seq)
        self.seqio = seqio
        self.parse_position()
        Diversity_Window.window_ids.append(self.id)
        Diversity_Window.lookup[self.name] = self


    def add_targets(self, startForward, endReverse):
        """add targets to object and create regex"""
        self.startForward = startForward
        self.endReverse = endReverse

        self.re_startForward, self.re_endReverse = \
            ['[ACTG]?'+target+'[ACTG]+' for target in [self.startForward, self.endReverse]]

    def matchToTarget(self, grepSeq):
        """match regex to grepSeq and return T/F"""
        if re.match(self.re_startForward, grepSeq):
            self.leftCuts += 1
            return True
        elif re.match(self.re_endReverse, grepSeq):
            self.rightCuts += 1
            return True
        return False

    def increment_highFidelityCount(self):
        """increment highFidelityCount"""
        self.highFidelityCount += 1

    def get_highFidelityCount(self):
        """class method to return the highFidelityCount"""
        return str(self.highFidelityCount)

    def get_totalCuts(self):
        """class method to calculate the total cuts"""
        return str(self.leftCuts + self.rightCuts)

    def get_leftCuts(self):
        """class method to return left-side cuts"""
        return str(self.leftCuts)

    def get_rightCuts(self):
        """class method to return right-side cuts"""
        return str(self.rightCuts)

    def parse_position(self):
        if self.id != None:
            dw, self.chrom, scope = self.id.split('::')
            self.left_pos, self.right_pos = [int(s) for s in scope.split('-')]
        self.name = self.id.split('::')[0]
        self.parse_chrom()
    def parse_chrom(self):
        # split species and pull chrom number (removing 0 from number if below 10)
        self.chrom_num = self.chrom.split('_')[1]
        if self.chrom_num[0] == '0':
            self.chrom_num = self.chrom_num.replace('0', '')
    def add_snp(self, snp):
        """add snp if the snp position is within the window scope"""
        snp_pos = int(snp.pos)
        if (snp_pos >= self.left_pos) and (snp_pos <= self.right_pos):
            self.snps.append(snp)
            snp.relative_position(self) # make relative position of snp on sequence
    def add_read(self, read):
        """check if read is unique and add to count of total and unique"""
        try:
            r_seq = str(read.seq)
            if r_seq not in self.unique_reads:
                self.unique_reads[r_seq] = 0
            self.unique_reads[r_seq] += 1
            self.total_reads += 1
        except AttributeError:
            # most likely from pairing failure
            pass
    def add_uniq(self, uniq):
        """add uniq dictionary to object and parse for total"""
        self.unique_reads = uniq
        for s, c in uniq.items():
            self.total_reads += c
    def add_mip(self, mip):
        """add a mip object to list of mips"""
        self.mips.append(mip)
    def add_snp_string(self, snp_str):
        """add snp string from alignment to list"""
        self.snp_strings.append(snp_str)
    def add_haplotype(self, haplo):
        """add haplotype snp string to unique dictionary and to list"""
        seq = str(haplo.seq)
        count = int(haplo.id.split('_')[1])
        if seq not in self.unique_haplotypes:
            self.unique_haplotypes[seq] = 0
        self.unique_haplotypes[seq] += count
        self.total_haplo += count
    def add_alignment_order(self, aln_order):
        """adds order of alignment and parses weights"""
        self.alignment_order = aln_order
        self.aln_order_weights = [int(s.split('_')[1]) for s in self.alignment_order]
    def snp_printout(self):
        """write to stdout"""
        for s in self.snps:
            info = [s.id, s.ref, s.nref, self.name]
            print('::'.join(info))
    def snp_write(self, output):
        """write to filename"""
        for s in self.snps:
            info = [s.id, s.ref, s.nref, self.name]
            output.write('::'.join(info) + '\n')
    def pass_threshold(self, seq, threshold):
        """True if the frequency of read is above threshold else False"""
        if (self.unique_reads[seq] / float(self.total_reads)) > threshold:
            return True
        return False
    def print_hf(self, threshold=0.001):
        """
        given a threshold as decimal of total reads print all reads with frequency
        above that threshold (0.001 == 0.1%)
        """
        print('>' + self.id)
        print(self.seq)
        for seq in self.unique_reads:
            if self.pass_threshold(seq, threshold) == True:
                print('>' + self.name + '_' + str(self.unique_reads[seq]))
                print(seq)
    def hf_gen(self, reference=True, threshold=0.001):
        """high fidelity read generator with reference sequence as first record"""
        reference = SeqRecord(Seq(self.seq), id=self.name)
        yield reference
        for seq in self.unique_reads:
            if self.pass_threshold(seq, threshold) == True:
                r = SeqRecord(Seq(seq), id=self.name + '_' + str(self.unique_reads[seq]))
                yield r
    def count_snps(self):
        """pass snp string to snp object to store information"""
        for i in range(len(self.snp_strings)):
            self.snps[i].count_snps(self.snp_strings[i], self.aln_order_weights)
    def haplotype_percentage(self, count):
        """return the percentage of a haplotype in total count"""
        return (float(count) / self.total_haplo * 100)
    def print_snp_count(self, ofn):
        """print the ratio of total, ref, and nref for each snp in csv format"""
        for s in self.snps:
            if s.get_distribution() != False:
                l = [s.full_position]
                [l.append(str(r)) for r in s.get_distribution()]
                ofn.write(','.join(l) + '\n')
    def print_haplotypes(self, ofn):
        """print the haplotype percentages for the window"""
        # total_haplotypes , percentage, haplotype_string, window
        for seq,count in self.unique_haplotypes.items():
            l = [self.total_haplo ,self.haplotype_percentage(count), seq, self.name]
            l_str = [str(i) for i in l]
            ofn.write(','.join(l_str) + '\n')
class SNP:
    snps = list()
    def __init__(self):
        self.chrom = None
        self.pos = None
        self.full_position = None
        self.id = None
        self.ref = None
        self.nref = None
        self.r_pos = None
        self.snp_string = None
        self.weights = None
        self.num_ref = 0
        self.num_nref = 0
        self.num_unknown = 0
        self.num_total = 0
        SNP.snps.append(self)
    def add_line(self, line):
        """pull attributes of each line and add aatributes to SNP object"""
        attrib = line.split('\t')
        self.chrom = str(int(attrib[0].split('_')[1]))
        self.pos, self.ref, self.nref = attrib[-4:-1]
        self.full_position = ':'.join([self.chrom, self.pos])

        # indexes = [0, 1, -6, -5, -4]
        # attrib = [line[i] for i in indexes]
        # self.chrom, self.pos, self.id, self.ref, self.nref = attrib
        # self.full_position = ':'.join([self.chrom, self.pos])
    def relative_position(self, dw):
        """make relative position of snp on reference sequence"""
        self.r_pos = int(self.pos) - int(dw.left_pos) - 1
    def count_snps(self, snp_string, weights):
        """gather relevant snp information and store within object"""
        self.snp_string = snp_string
        self.weights = weights
        for i in range(len(snp_string)):
            c = snp_string[i]
            if c == self.ref:
                self.num_ref += weights[i]
            elif c == self.nref:
                self.num_nref += weights[i]
            elif c != '-':
                self.num_unknown += weights[i]
            if c != '-':
                self.num_total += weights[i]
    def percentage(self, num):
        return float(num)/self.num_total * 100
    def get_distribution(self):
        """return: total, percent_ref, percent_nref, unknown"""
        if self.num_total != 0:
            return self.num_total, self.percentage(self.num_ref), self.percentage(self.num_nref), self.percentage(self.num_unknown)
        return False
class Alignment:
    alignments = list()
    def __init__(self):
        self.alignment = None
        self.reference = None
        self.order_of_id = list()
        self.snp_locations = list()
        self.snps = list()
        Alignment.alignments.append(self)
    def mips_and_dw(self, dw):
        """align mips to dw: also add snps to Alignment class for easy usage"""
        i_fn, o_fn = self.make_alignment_fasta(dw)
        self.align(i_fn, o_fn)
        [self.snp_locations.append(s.r_pos) for s in dw.snps]
        [self.snps.append(s) for s in dw.snps]
    def hf_and_dw(self, gen, dw):
        """align high fidelity reads to dw, also adds snps to Alignment class"""
        self.alignment = self.generator_alignment(gen)
        self.reference = self.alignment[0]
        dw.add_alignment_order([s.id for s in self.alignment[1:]])
        [self.snp_locations.append(s.r_pos) for s in dw.snps]
        [self.snps.append(s) for s in dw.snps]
    def generator_alignment(self, gen):
        """perform clustalo alignment given a generator (will make temp file and return AlignIO object)"""
        # run alignment on generated reads
        proc = subprocess.Popen(['./clustalo --outfmt=clu --infile=-'], shell=True, stdin=PIPE, stdout=PIPE)

        # write to process stdin
        [proc.stdin.write('>' + r.id + '\n' + str(r.seq) + '\n') for r in gen]


        # store stdout
        stdout, stderr = proc.communicate()

        # create temporary file to store stdout to pass to alignIO
        temp = tempfile.TemporaryFile()
        temp.write(stdout)
        temp.seek(0)

        # return AlignIO object
        return AlignIO.read(temp, 'clustal')
    def make_alignment_fasta(self, dw):
        """create a fasta containing diversity window and mips associated"""
        i_fn = dw.name + '_alignment_query.fasta'
        o_fn = dw.name + '_alignment_output.fasta'
        fasta = Fasta()
        fasta.dw_write(i_fn, dw)
        [fasta.mip_write(i_fn, m, direction=True) for m in dw.mips]
        return i_fn, o_fn
    def align(self, i_fn, o_fn):
        """call clustalo and align query"""
        subprocess.call(['clustalo -i ' + i_fn + ' -o ' + o_fn + ' --outfmt=clu --force'], shell=True)
        self.alignment = AlignIO.read(open(o_fn, 'r'), 'clustal')
        self.reference = self.alignment[0]
        self.remove(i_fn)
        self.remove(o_fn)
    def remove(self, fn):
        """delete a filename from directory"""
        subprocess.call(['rm ' + fn], shell=True)
    def get_snps_columns(self):
        """return expected genotype of snp positions excluding reference"""
        self.update_snps()
        return [self.alignment[:,i][1:] for i in self.snp_locations]
    def get_snps_rows(self):
        """return only the snp strings for each unique read"""
        for a in self.alignment[1:]:
            seq = ''.join([a.seq[i] for i in self.snp_locations])
            r = SeqRecord(Seq(seq), id=a.id)
            yield r
    def gap_locations(self):
        """check the diversity window for gaps and append to gap_locations"""
        return [i for i in range(len(self.reference.seq)) if self.reference.seq[i] == '-']
    def update_snps(self):
        """update relative position of snps by gap locations"""
        for g in self.gap_locations():
            for s in self.snp_locations:
                if g <= s:
                    index = self.snp_locations.index(s)
                    # increase remaining reference positions in snps by one
                    for i in range(index, len(self.snp_locations)):
                        self.snp_locations[i] += 1
                    break # necessary break
class Pair:
    pairs = list()
    lookup = dict()
    uniq = dict()
    uniq_total = 0
    def __init__(self, uid):
        self.uid = uid
        self.left = None
        self.right = None
        self.count = 0
        Pair.lookup[uid] = self
        Pair.pairs.append(self)
    def add_seq(self, seqio):
        """add sequence to object and join when both pairs are added"""
        direction = seqio.id.split('_')[1]
        if direction == '1':
            self.left = seqio
        else:
            self.right = seqio
        self.count += 1
        if self.count == 2:
            self.join_pairs()
            return True
    def join_pairs(self, print_bool=False):
        """join the two sequences on the longest overlap"""
        s1 = str(self.left.seq)
        s2 = str(self.right.seq)
        ss = self.get_overlap(s1, s2)
        try:
            s1, s2, ss = self.find_direction(s1, s2, ss)
            self.join(s1,s2,ss)
            if print_bool == True:
                self.print_overlap(s1,s2,ss)
        except TypeError:
            pass
    def get_overlap(self, s1, s2):
      """return the largest region of overlap between two seqs"""
      s = difflib.SequenceMatcher(None, s1, s2)
      pos_a, pos_b, size = s.find_longest_match(0, len(s1), 0, len(s2))
      return s1[pos_a:pos_a+size]
    def find_direction(self, s1, s2, ss):
        """return strings in proper orientation"""
        if s1.index(ss) < s2.index(ss):
            return [s2, s1, ss]
        elif s1.index(ss) > s2.index(ss):
            return [s1, s2, ss]
    def print_overlap(self, s1, s2, ss):
        """format print overlap"""
        print(s1)
        print('.' * s1.index(ss) + ss)
        print('.' * s1.index(ss) + s2)
    def join(self, s1, s2, overlap):
        """s1 until overlap, overlap, s2 past overlap"""
        index1 = s1.index(overlap) # position of overlap in s1
        index2 = s2.index(overlap) # position of overlap in s2
        return (s1[:index1] + s2[index2:])
        self.append_uniq()
    def append_uniq(self):
        """append sequence to unique dictionary and add to total"""
        if self.seq not in Pair.uniq:
            Pair.uniq[self.seq] = 0
        Pair.uniq[self.seq] += 1
        Pair.uniq_total += 1
