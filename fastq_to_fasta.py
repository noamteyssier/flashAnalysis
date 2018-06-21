import Bio, sys, os, argparse, argcomplete
from Bio import SeqIO

def get_filenames(directory):
    """return the filenames of the fast(q/a)s or complain directory doesnt exist"""
    ext = '.fq'
    try:
        g = [directory+f for f in os.listdir(directory) if ext in f]
        if len(g) == 2:
            return g
        else:
            print 'Must only be 2 <*%s> in directory <%s>' %(ext, directory)
    except OSError:
        print 'Cannot Find Given Directory <%s>' %directory
    sys.exit()
def q2a(filename, direction, f):
    with open(filename, 'r') as fastq:
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
    with open(output, 'w+') as f:
        if '.fasta' not in forward:
            q2a(forward, 'F', f)
        if '.fasta' not in reverse:
            q2a(reverse, 'R', f)
def __main__():
    parser = argparse.ArgumentParser(description='convert two illumina fastq files into three fastas (forward, reverse, and both)')
    parser.add_argument('-d', '--directory', help='directory containing fastq files (can only have two)', required=True)
    parser.add_argument('-o', '--output', help='output filename for converted file', required=True)
    argcomplete.autocomplete(parser)
    args = parser.parse_args()

    # assign output filename
    output_fn = args.output if '.fasta' in args.output else args.output+'.fasta'

    # perform fastq -> fasta
    fn1, fn2 = get_filenames(args.directory)
    write_fastas(fn1, fn2, output_fn)


__main__()
