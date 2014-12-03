"""
Generates a numpy matrix of labels given an
intersection of bed sequences and genomes.
Serializes numpy array to outfile.
"""

import sys
from Bio import SeqIO
from features import *
import numpy as np


def hack_key(key):
    chrom, start, end = key.strip().split("\t")
    start = start[:-2]
    end = end[:-2]
    return "\t".join([chrom, start, end]) + "\n"


def get_key(line):
    key = "\t".join(line.split()[:3]) + "\n"
    return hack_key(key)


def get_val(line):
    return line.split()[6]


def gen_tf_dict(fin):
    # can't open as regular bedfile; instead, have to intersect
    f = open(fin)
    uniq_keys = set([get_key(line) for line in f])

    result = dict([(k, []) for k in uniq_keys])
    f.seek(0, 0)
    for line in f:
        key = get_key(line)
        assert key in uniq_keys
        result[key].append(get_val(line))
    return result




if __name__ == "__main__":
    if len(sys.argv) < 5:
        print """
        Arguments: <intersection.bed> <examples.fa> <out.npy> <op>
            First is intersection bed file with region info from both original bed files (-wa -wb)
            Second is path to bedfile
            Third is outfile.
            fourth is the operation. Can be (tf|<something else>)
        """
        sys.exit(1)

    fintersect = sys.argv[1]
    fexamples = sys.argv[2]
    fout = sys.argv[3]
    op = sys.argv[4]

    tfdict = gen_tf_dict(fintersect)
    seqfile = [x for x in SeqIO.parse(fexamples, "fasta")]
    labels = []

    for x in seqfile:
        if op == "tf":
            label = np.zeros(3, dtype=np.int32)
            key = hack_key(description_to_bed(x.description))
            if key in tfdict:
                label = extract_feat_tf(tfdict[key])
        else:
            # annotated fasta file, so can use seq_to_bed
            label = np.zeros(4, dtype=np.int32)
            if op in x.description:
                key = hack_key(seq_to_bed(x.description))
                if key in tfdict:
                    label = extract_feat_nontf(tfdict[key])
        labels.append(label)

    result = np.array(labels)
    np.save(fout, result)
