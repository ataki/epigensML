""" 
Nice plots
"""
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve
from sklearn.preprocessing import label_binarize
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt


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
    diff = np.array([1 if y_true[i] == preds[i] else 0 for i in range(len(y_true))])
    cs = "rg"
    cats = [0, 1]
    target_names = ["incorrect", "correct"]
    for c, i, target_name in zip(cs, cats, target_names):
        plt.scatter(X_r[diff == i, 0], X_r[diff == i, 1], c=c, label=target_name)
        plt.legend()
        plt.title("PCA of correct/incorrect predictions")
    # plt.show()
    plt.savefig("figures/residual-scatter.png")


def plot_precision_recall(y_true, y_scores):
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)
    plt.figure()
    plt.plot(recall, precision, 'g-')
    plt.title("Precision-Recall Curve")
    plt.savefig("figures/pr-curve.png")


def plot_roc(y_true, y_score, title, out="figures/roc-curve.png"):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    is_1d = len(y_true.shape) == 1

    # Plot a curve for each label
    if is_1d: 
        fpr[0], tpr[0], _ = roc_curve(y_true, y_score)
        roc_auc[0] = auc(fpr[0], tpr[0])
    else:
        for i in range(y_true.shape[1]):
            fpr[i], tpr[i], _ = roc_curve(y_true[:, i], y_score[:, i])
            roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area
    fpr["micro"], tpr["micro"], _ = roc_curve(y_true.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    # Plot of a ROC curve for a specific class
    for i in fpr.keys(): 
        if i != "micro":
            plt.figure()
            plt.plot(fpr[i], tpr[i], label='ROC curve (area = %0.2f)' % roc_auc[i])
            plt.plot([0, 1], [0, 1], 'k--')
            plt.xlim([0.0, 1.0])
            plt.ylim([0.0, 1.05])
            plt.xlabel('False Positive Rate')
            plt.ylabel('True Positive Rate')
            plt.title(title)

    plt.legend(loc="lower right")
    # plt.show()
    plt.savefig(out)

