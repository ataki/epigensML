#!/usr/bin/python

"""
Given two files "pos.fa" and "neg.fa", does 5-fold
cross validation to determine the accuracy of predicting
positive from negative examples.
Optionally writes the model to an "out" directory.
Optionally plots precision/recall curves.
"""

import re
import sys
import pdb
import itertools
import pandas as pd
import numpy as np
import pybedtools

from Bio import SeqIO
from sklearn import svm
from sklearn import cross_validation
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.externals import joblib
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold

from matplotlib import pyplot as plt

# === Config ===

FEATURE_SELECTION = False

FOLD_CV = False

SPLIT_CV = True

NORMALIZE = True

GLOBAL_K = 6


# === Preprocessing ===


def reverse(input_str):
    """ Simple reverse string """
    list_form = list(input_str)
    list_form.reverse()
    return list_form


def cmpl_base_pair(x):
    """ Get complementary base pair """
    if x == 'A': return 'T'
    elif x == 'C': return 'G'
    elif x == 'T': return 'A'
    elif x == 'G': return 'C'
    else: return 'N'


def neg_strand(pos_strand):
    """ Given pos_strand sequence, returns complementary
    negative sequence """
    return "".join(map(cmpl_base_pair, reverse(pos_strand)))


def gen_kmers(seq, k=GLOBAL_K):
    """ Returns generator for all kmers in given DNA sequence.
    May contain duplicates """
    return (seq[i:i+k] for i in range(len(seq) - k + 1))


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


def get_locations_to_y_tIndex(locations):
    """ mocations_to_y_tIndex is a dictionary that maps location
    (eg. hindbrain) to  indices into the y_t vector. """

    locations_to_y_tIndex = {
        'forebrain': [],
        'hindbrain': [],
        'limb': [],
        'rest': []
    }

    cutoff = (8 * len(locations)) / 10

    index = 0
    for (x, y) in locations[cutoff:]:
        if len(y) > 0:
            for location in y:
                if "forebrain" in location:
                    locations_to_y_tIndex['forebrain'].append(index)
                if "hindbrain" in location:
                    locations_to_y_tIndex['hindbrain'].append(index)
                if "limb" in location:
                    locations_to_y_tIndex['limb'].append(index)
                if "forebrain" not in location and \
                        "hindbrain" not in location and \
                        "limb" not in location:
                    if index not in locations_to_y_tIndex['rest']:
                        locations_to_y_tIndex['rest'].append(index)
        index += 1
    return locations_to_y_tIndex


"""
Dist:
    forebrain => 309
    midbrain => 253
    hindbrain => 236
    neural tube => 170
    limb => 162
    other => 596
"""


def iftl(tissue_label):
    """ Returns int for tissue label """
    if tissue_label == "brain": return 1
    elif tissue_label == "limb": return 2
    elif tissue_label == "neural": return 3
    else: return 0


def lftd(description):
    """ converts descriptions to tissue labels.
    Input is fasta descr. """
    regex = "([^\[]+)\[(\d+)\/(\d+)\]"
    ranks = []
    for raw in description.split("|")[4:]:
        line = raw.strip()
        if "brain" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("brain", int(count)))
        if "limb" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("limb", int(count)))
        if "neural" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("neural", int(count)))
    if len(ranks) == 0:
        return 0
    else:
        label, score = max(ranks, key=lambda x: x[1])
        return iftl(label)


def ifbd(brain_label):
    if brain_label == "fore": return 1
    elif brain_label == "mid": return 2
    elif brain_label == "hind": return 3
    return 0


def lfbd(description):
    """ converts descriptions to brain labels.
    Input is fasta descr """
    regex = "([^\[]+)\[(\d+)\/(\d+)\]"
    ranks = []
    for raw in descriptions.split("|")[4:]:
        line = raw.strip()
        if "midbrain" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("mid", int(count)))
        if "forebrain" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("fore", int(count)))
        if "hindbrain" in line:
            _, count, _ = re.match(regex, line).groups()
            ranks.append(("hind", int(count)))
    if len(ranks) == 0:
        return 0
    else:
        label, score = max(ranks, key=lambda x: x[1])
        return ifbd(label)


def parse_fa(path, label):
    """ Given a fasta file that represents a label
    class, returns a pair of (sequence, label) numpy
    arrays. Useful for constructing X/y for training """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')
    seqs = []
    _labels = []

    for entry in human_fasta_seq:
        seqs.append(str(entry.seq).lower())
        _labels.append(float(label))

    return (seqs, _labels)


def parse_fa_tissue(path):
    """ Given a fasta file that represents positive
    labels, returns a pair of (sequence, label) numpy
    arrays. Useful for constructing X/y for training """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')

    seqs = []
    _labels = []

    for entry in human_fasta_seq:
        seqs.append(str(entry.seq).lower())
        _labels.append(lftd(entry))

    return (seqs, _labels)


def parse_fa_fine_grain(path):
    """ Given a fasta file that represents positive
    labels, returns a pair of (sequence, label) numpy
    arrays. Useful for constructing X/y for training """
    fasta_file = open(path)
    human_fasta_seq = SeqIO.parse(fasta_file, 'fasta')
    seqs = []
    _labels = []

    for entry in human_fasta_seq:
        seqs.append(str(entry.seq).lower())
        _labels.append(lfbd(entry))

    return (seqs, _labels)


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
    y = np.array(labels)
    print "train. matrix dims (X): ", X.shape
    print "num labels (y): ", len(y)
    print "+ ", len(np.where(y == 1)[0])
    print "- ", len(np.where(y == -1)[0])
    print "------------------------------------"
    return (X, y)


def get_TAAT_core_col(examples):
    all_seqs = [x[1] for x in examples]
    regex = "^([atgc])+(taat)([atcg])+$"
    expr = lambda x: 1.0 if re.match(regex, x) else 0.0
    return np.array(map(expr, all_seqs)).reshape(len(all_seqs), 1)


def get_Ebox_col(examples):
    all_seqs = [x[1] for x in examples]
    regex = "ca[atcg]{2}tg"
    expr = lambda x: 1.0 if re.search(regex, x) else 0.0
    return np.array(map(expr, all_seqs)).reshape(len(all_seqs), 1)


# === Prediction ===


def plot_2d_results(X, y, preds):
    pca = PCA(n_components=2)
    X_r = pca.fit(X).transform(X)

    # Plot scatter
    plt.figure()
    cs = "cm"
    cats = [1, -1]
    target_names = ["positive", "negative"]
    for c, i, target_name in zip(cs, cats, target_names):
        plt.scatter(X_r[y == i, 0], X_r[y == i, 1], c=c, label=target_name)
    plt.legend()
    plt.title("PCA of 2d data")
    plt.savefig("figures/data-scatter.png")

    # Plot mispredictions
    plt.figure()
    diff = np.array([1 if y_test[i] == preds[i] else 0 for i in range(len(y_test))])
    cs = "rg"
    cats = [0, 1]
    target_names = ["incorrect", "correct"]
    for c, i, target_name in zip(cs, cats, target_names):
        plt.scatter(X_r[diff == i, 0], X_r[diff == i, 1], c=c, label=target_name)
        plt.legend()
        plt.title("PCA of correct/incorrect predictions")
    # plt.show()
    plt.savefig("figures/residual-scatter.png")


def plot_precision_recall(y_test, y_scores):
    precision, recall, thresholds = precision_recall_curve(y_test, y_scores)
    plt.figure()
    plt.plot(recall, precision, 'g-')
    plt.title("Precision-Recall Curve")
    plt.savefig("figures/pr-curve.png")


def plot_roc(y_test, y_score):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(1):
        fpr[i], tpr[i], _ = roc_curve(y_test, y_score)
        roc_auc[i] = auc(fpr[i], tpr[i])
    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    # Plot of a ROC curve for a specific class
    plt.figure()
    plt.plot(fpr[0], tpr[0], label='ROC curve (area = %0.2f)' % roc_auc[0])
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC, Kmer counts used to predict general enhancer functionality')
    plt.legend(loc="lower right")
    # plt.show()
    plt.savefig("roc-curve.png")


def print_usage_and_exit():
    print "Usage: enhancer_clf.py pos.fa neg.fa <prediction_type>"
    print "First two args are paths to datasets"
    print "[Optional] third arg is one of <enhancer|tissue|fine-grain>"
    print "enhancer predicts general enhancer activity"
    print "tissue predicts limb/brain/heart"
    print "fine-grain predicts forebrain/hindbrain"


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print_usage_and_exit()

    pos_dataset = sys.argv[1]
    neg_dataset = sys.argv[2]
    prediction_type = sys.argv[3] if len(sys.argv) > 3 else "enhancer"

    if prediction_type == "enhancer":
        pos_seq, pos_labels = parse_fa(pos_dataset, 1)
        neg_seq, neg_labels = parse_fa(neg_dataset, -1)
        examples, labels = pos_seq + neg_seq
        labels = pos_labels + neg_labels
    elif prediction_type == "tissue":
        examples, labels = parse_fa_tissue(pos_dataset)
    elif prediction_type == "fine-grain":
        examples, labels = parse_fa_fine_grain(pos_dataset)

    # feature vector index :=> kmer string
    kmers_index = get_kmers_index_lookup()

    # feature matrix and label vector
    X, y = get_XY(examples, labels, kmers_index)

    # scale raw counts
    if NORMALIZE:
        X = normalize(X, axis=1, norm='l1')

    # Add e-box and taat core cols
    ebox_col = get_Ebox_col(examples)
    taat_col = get_TAAT_core_col(examples)
    X = np.hstack((X,ebox_col, taat_col))

    clf = svm.SVC(kernel='rbf', C=1)

    if FEATURE_SELECTION:
        print "Feature selecting top 10 features"
        from sklearn.feature_selection import SelectKBest
        from sklearn.feature_selection import chi2
        # Remove low-variance features
        # K-best features
        X = SelectKBest(chi2, k=10).fit_transform(X, y)

    if FOLD_CV:
        print "Performing 5-fold cv"
        scores = cross_validation.cross_val_score(clf, X_train, y_train, cv=5)
        print "%d-fold cv, average accuracy %f" % (len(scores), scores.mean())

    if SPLIT_CV:
        print "Performing train/test split cv"
        X_train, X_test, y_train, y_test = cross_validation.train_test_split(
            X, y, test_size=0.3, random_state=0)
        clf.fit(X_train, y_train)
        clf.score(X_test, y_test)

        # transform labels from [-1,1] to [0,1]
        _y_test = label_binarize(y_test, classes=[-1, 1])
        y_scores = clf.decision_function(X_test)

        print "Plotting results"
        plot_roc(_y_test, y_scores)
        plot_precision_recall(_y_test, y_scores)
        plot_2d_results(X_test, y_test, clf.predict(X_test))
