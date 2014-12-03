"""
Base script. Intended to be forked and modified.

Given two files "pos.fa" and "neg.fa", does 5-fold
cross validation to determine the auroc predicting
positive from negative examples.

Optionally writes the model to an "out" directory.

Optionally plots precision/recall curves.
"""

# Imports for IPython
import pdb
import itertools
import argparse
import collections
import numpy as np
from Bio import SeqIO

# make util functions available
from features import *
from plot import *

# Script imports
from sklearn import svm
from sklearn import ensemble
from sklearn import cross_validation as cv
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from sklearn.multiclass import OneVsRestClassifier
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2

from matplotlib import pyplot as plt

FEATURE_SELECTION = False

FOLD_CV = False

PLOT_RESULTS = not FOLD_CV

EXTRA = True

NORMALIZE = True

import time
if __name__ == "__main__":
    start = time.clock()
    parser = argparse.ArgumentParser()
    parser.add_argument("pos_exs", help="path to pos examples")
    parser.add_argument("neg_exs", help="path to neg examples")
    parser.add_argument("--pos_tf", help="path to preprocessed intersection labels npy file")
    parser.add_argument("--neg_tf", help="path to preprocessed intersection labels npy file")
    args = parser.parse_args()

    pos_dataset = args.pos_exs
    neg_dataset = args.neg_exs
    pos_tf = args.pos_tf if args.pos_tf else None
    neg_tf = args.neg_tf if args.neg_tf else None

    pos_seq, pos_labels = parse_fa(pos_dataset, 1)
    neg_seq, neg_labels = parse_fa(neg_dataset, -1)
    pos_extra = np.load(pos_tf) if pos_tf else None
    neg_extra = np.load(neg_tf) if neg_tf else None

    # Vertically stack everything
    examples = np.concatenate((pos_seq, neg_seq))
    labels = np.concatenate((pos_labels, neg_labels))
    extra_feats = np.concatenate((pos_extra, neg_extra)) if pos_extra is not None and neg_extra is not None else None
    extra_feats = extra_feats[:, 0:2]

    # feature vector index :=> kmer string
    kmers_index = get_kmers_index_lookup()

    # feature matrix and label vector
    X, y = get_XY(examples, labels, kmers_index)

    # scale raw counts
    if NORMALIZE:
        X = normalize(X, axis=1, norm='l1')

    # == Add extra features ==

    # Add Kristin's features
    if EXTRA:
        if extra_feats is not None:
            X = np.hstack((X, extra_feats))

    # Add e-box and taat core cols
    # ebox_col = get_ebox_col(examples)
    # taat_col = get_taat_col(examples)
    # X = np.hstack((X, taat_col))

    clf = svm.SVC(kernel='linear')

    if FEATURE_SELECTION:
        print "Feature selecting top 10 features"
        # Remove low-variance features
        # K-best features
        X = SelectKBest(chi2, k=600).fit_transform(X, y)

    if FOLD_CV:
        print "Performing 5-fold cv"
        scores = cv.cross_val_score(
            clf, X, y, cv=5, scoring="roc_auc"
        )
        print "%d-fold cv, average auRoc %f" % (len(scores), scores.mean())

    if PLOT_RESULTS:
        print "Performing train/test split cv"
        X_train, X_test, y_train, y_test = cv.train_test_split(
            X, y, test_size=0.3, random_state=0
        )
        clf.fit(X_train, y_train)

        print "Plotting results"
        y_scores = clf.decision_function(X_test)
        if EXTRA:
            plot_roc(y_test, y_scores, "ROC Enhancer",
                out="figures/roc-curve-enh-fsel-extra.png")
        else:
             plot_roc(y_test, y_scores, "ROC Enhancer",
                out="figures/roc-curve-combo-vs-rand.png")
       # plot_precision_recall(y_true, y_scores)
        # plot_2d_results(X_test, y_test, clf.predict(X_test))
        print "Done plotting"

    end = time.clock()
    print "time taken : %fs" % (end - start)
