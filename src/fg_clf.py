"""
Given an annotated fasta file of brain enhancer sequences,
trains and evaluates a classifier to distinguish among
different locations of the brain.

Classifiers used are one-vs-one to distinguish btwn
two tissues, and one-vs-all for the rest.
"""

# Imports for IPython
import argparse
import numpy as np

# make util functions available
from features import *
from plot import *

# Script imports
from sklearn import svm
from sklearn import cross_validation as cv
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import normalize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

FEATURE_SELECTION = False

FOLD_CV = False

PLOT_RESULTS = not FOLD_CV

NORMALIZE = True

GLOBAL_K = 6


def filter_one_v_all(description):
    """ Filters descriptions for those that have
    the necessary parts """
    brain_parts = ["forebrain", "midbrain", "hindbrain"]
    for part in brain_parts:
        if part in description:
            return True
    return False


def filter_binary(description):
    return "brain" in description


def one_v_one(description):
    """ Treats the first two entries in brain-parts
    as the labels for a one-v-one clf """
    return "forebrain" in description


def one_v_all(description):
    """ converts descriptions to multi-label
    brain label.
    Indices in label correspond to:
        1 => "forebrain"
        2 => "midbrain"
        3 => "hindbrain"
    """
    label = [0, 0, 0]
    for raw in description.split("|")[4:]:
        line = raw.strip()
        if "midbrain" in line:
            label[0] = 1
        if "forebrain" in line:
            label[1] = 1
        if "hindbrain" in line:
            label[2] = 1
    return label


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pos_exs", help="path to pos examples")
    parser.add_argument("--forebrain", help="path to forebrain npy file")
    parser.add_argument("--prediction_type", help="binary or multilabel", default="multilabel")
    args = parser.parse_args()

    pos_dataset = args.pos_exs
    prediction_type = args.prediction_type
    forebrain_feats = None

    if prediction_type == "multilabel":
        examples, labels, kept_rows = parse_fa_fine_grain(
            pos_dataset, one_v_all, filter_one_v_all
        )
    else:
        examples, labels, kept_rows = parse_fa_fine_grain(
            pos_dataset, one_v_one, filter_binary
        )

    if args.forebrain:
        forebrain_feats = np.load(args.forebrain)

    # feature vector index :=> kmer string
    kmers_index = get_kmers_index_lookup()

    # feature matrix and label vector
    X, y = get_XY(examples, labels, kmers_index)

    # scale raw counts
    if NORMALIZE:
        X = normalize(X, axis=1, norm='l1')

    # Add e-box and taat core cols
    # ebox_col = get_ebox_col(examples)
    # taat_col = get_taat_col(examples)
    # X = np.hstack((X, taat_col))

    # == Extra features ==
    if forebrain_feats is not None:
        feats = forebrain_feats[kept_rows, :]
        X = np.hstack((X, feats))

    clf = None

    if prediction_type == "multilabel":
        svc = svm.SVC(kernel='rbf')
        clf = OneVsRestClassifier(svc)
    else:
        clf = svm.SVC(kernel='rbf')

    if FEATURE_SELECTION:
        print "Feature selecting top 10 features"
        # Remove low-variance features
        # K-best features
        X = SelectKBest(chi2, k=10).fit_transform(X, y)

    print "Performing 5-fold cv"
    scores = cv.cross_val_score(
        clf, X, y, cv=5, scoring="roc_auc"
    )
    print "%d-fold cv, average auRoc %f" % (len(scores), scores.mean())
