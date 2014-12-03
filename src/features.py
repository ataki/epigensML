"""
Feat. extraction for parts 1/2
"""
import re
import pdb
import numpy as np
import itertools
import pybedtools
from Bio import SeqIO

GLOBAL_K = 6

brain_seqFile = None
heart_seqFile = None
limb_seqFile = None

# === Part 1: Kmer classifiers ===


def reverse(input_str):
    """ Simple reverse string """
    list_form = list(input_str)
    list_form.reverse()
    return list_form


def cmpl_base_pair(x):
    """ Get complementary base pair """
    if x == 'A':
        return 'T'
    elif x == 'C':
        return 'G'
    elif x == 'T':
        return 'A'
    elif x == 'G':
        return 'C'
    else:
        return 'N'


def neg_strand(pos_strand):
    """ Given pos_strand sequence, returns complementary
    negative sequence """
    return "".join(map(cmpl_base_pair, reverse(pos_strand)))


def gen_kmers(seq, k=GLOBAL_K):
    """ Returns generator for all kmers in given DNA sequence.
    May contain duplicates """
    return (seq[i:i + k] for i in range(len(seq) - k + 1))


def locfd(description):
    """ converts descriptions to locations. Input
    is a fasta descr. delimited by "|". Whether
    the region is enhanced is denoted by the ?th
    entry in the array """
    return description.split("|")[4:]


def lts(label):
    """ converts label to sign """
    return 1 if label == "positive" else -1


def lfd(description):
    """ converts descriptions to labels. Input
    is a fasta descr. delimited by "|". Whether
    the region is enhanced is denoted by the 4th
    entry in the array """
    return lts(description.split("|")[3].strip())


def description_to_bed(description):
    return description.split()[1].strip() \
        .replace("range=", "") \
        .replace(":", "\t") \
        .replace("-", "\t") + "\n"


def parse_fa(path, label):
    """ Given a fasta file that represents a label
    class, returns sequence, label, and aux numpy
    arrays. Useful for constructing X/y for training """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')
    seqs = []
    labels = []

    for entry in human_fasta_seq:
        seqs.append(str(entry.seq).replace("n", "").lower())
        labels.append(float(label))

    seqs = np.array(seqs)
    labels = np.array(labels)
    return (seqs, labels)


def parse_fa_tissue(path, label_fn, filter_fn):
    """
    Given
        fasta file that represents positive labels
        label extracting fn given description

    Returns
        a pair of (sequence, label) numpy
        arrays. Useful for constructing X/y for training
    """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')

    seqs = []
    labels = []
    kept_rows = []

    idx = 0
    for entry in human_fasta_seq:
        if filter_fn(entry.description):
            seqs.append(str(entry.seq).replace("n", "").lower())
            labels.append(label_fn(entry.description))
            kept_rows.append(idx)
        idx += 1

    seqs = np.array(seqs)
    labels = np.array(labels)
    kept_rows = np.array(kept_rows)
    return (seqs, labels, kept_rows)


def parse_fa_fine_grain(path, label_fn, filter_fn):
    """
    Given
        fasta file that represents positive labels
        label extracting fn given description

    Returns
        a pair of (sequence, label) numpy
        arrays. Useful for constructing X/y for training
    """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')
    seqs = []
    labels = []
    kept_rows = []

    idx = 0
    for entry in human_fasta_seq:
        if filter_fn(entry.description):
            seqs.append(str(entry.seq).replace("n", "").lower())
            labels.append(label_fn(entry.description))
            kept_rows.append(idx)
        idx += 1

    seqs = np.array(seqs)
    labels = np.array(labels)
    kept_rows = np.array(kept_rows)
    return (seqs, labels, kept_rows)


def get_kmer_counts(seq, ref):
    """ Given example sequence and a reference table mapping
    kmers to indices, returns a numpy array representing one row
    of the feature vector.

    NOTE Finds kmers on both strands, avoids double-counting.
    This is an assumption made due to Figure 1 of the Lee et al. paper
    which describes a feature vector x with such counts:

        5' AAAAAA 3'  |>> x1
        3' TTTTTT 5'  |
        -------------------------
        ...
        -------------------------
        5' TTTAAA 3'  |>> xn
        3' AAATTT 5'  |
        -------------------------

    Based on the last entry, if "TTTAAA" is the same on both strands,
    don't count it.
    """
    row = np.zeros(len(ref))
    pos_kmers = gen_kmers(seq)
    for kmer in pos_kmers:
        if kmer != neg_strand(kmer):
            idx = ref[kmer]
            row[idx] += 1
    return row


def make_index_dict_from_list(l):
    """ Given list, creates dictionary where keys
    are contents of array and value is index of
    original array """
    return dict([(x, i) for i, x in enumerate(l)])


def get_kmers_index_lookup():
    """ Builds up mapping of index to a unique kmer """
    global GLOBAL_K
    all_kmers = [''.join(x) for x in itertools.product("atcg", repeat=GLOBAL_K)]
    return make_index_dict_from_list(list(set(all_kmers)))


def get_XY(examples, labels, kmer_index):
    X = np.vstack([get_kmer_counts(x, kmer_index) for x in examples])
    y = np.array(labels, dtype=np.int32)
    print "train. matrix dims (X): ", X.shape
    print "num labels (y): ", len(y)
    if len(labels.shape) == 1:
        print "+ ", len(np.where(y == 1)[0])
        print "- ", len(np.where(y != 1)[0])
    print "------------------------------------"
    return (X, y)


def get_taat_col(examples):
    regex = "^([atgc])+(taat)([atcg])+$"
    expr = lambda x: 1.0 if re.match(regex, x) else 0.0
    return np.array(map(expr, examples)).reshape(len(examples), 1)


def get_ebox_col(examples):
    regex = "ca[atcg]{2}tg"
    expr = lambda x: 1.0 if re.search(regex, x) else 0.0
    return np.array(map(expr, examples)).reshape(len(examples), 1)


# == Part 2: Additional features ==


def seq_to_bed(description):
    description = description.split("|")[1]
    chromName = description.split(":")[0]

    chromRange = description.split(":")[1]
    chromStart = chromRange.split("-")[0]
    chromEnd = chromRange.split("-")[1]

    bedLine = [chromName, chromStart, chromEnd]
    bedLine = '\t'.join(bedLine) + "\n"
    return bedLine


def extract_feat_tf(annotations):
    """ Given a bedtools line representing a sequence,
    returns row of general enhancer features """
    joint = "".join(annotations)
    features = [
        ('P300' in joint),
        ('TCF' in joint),
        ('TBF' in joint)
    ]
    return np.array(features, dtype=np.int32)


def extract_feat_nontf(annotations):
    """ Given single sequence that's labeled with
    a part of the brain, returns a list of additional
    epigenetic features """
    joint = "".join(annotations)
    features = [
        ('Enh' in joint),
        ('EnhG' in joint),
        ('Het' in joint),
        ('TxWk' in joint)
    ]
    return np.array(features, dtype=np.int32)
