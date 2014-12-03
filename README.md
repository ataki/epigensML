epigensML
=========

CS 273a: Predicting Tissue Specific Enhancer Activity from Epigenetic Marks and Sequence

## Authors

Nishith K., Kristin M., Jim Z. ([ataki](https://github.com/ataki))

## Summary

For reference, we're using the paper by Lee et al. in papers/2197 as a guide.
Their analysis uses only genome sequence features plus a set of general
sequence features in order to detect enhancer regions.

They use an SVM with hard boundary to separate +/- regions. Their SVM was
immune to kernel choice, which indicated linear separability across a 
high dimensional dataset (feature space of all possible kmers in regions of
interest).

## What's in Here

- `tasks` folder contains main prediction scripts.
- `src` folder contains python code for classifiers
- `scripts` folder contains data extraction code
- `requirements.txt` describes all our project dependencies, sans bedtools.

## Description of Data

Mostly for reference.

- Vista Dataset. Main one to use. Includes hg19/mm9 data, pos/neg enhancer,
  tissue, and part of brain labels.
- Beer Labs Dataset. Obtained >.90 clf accuracy for some pos/neg enhancer regions.
  Useful for testing that our models aren't too far off.

## Tasks

1. Predict general enhancer activity
2. Predict enhancer activity for tissue type.
3. Predict enhancer activity for parts of tissue.

## Techniques

We used an SVM because

- theoretical guarantees against overfitting
- historically proven model
- does well for small datasets

Our evaluation metric was au-ROC, which gives the probability that a randomly
chosen positive example ranks higher than a randomly chosen negative example.

We used one-vs-one to break up multi-class classification problems. For
multi-label prediction tasks, we used one-vs-rest.

In addition to kmer counts, we boosted Task 1 by adding indicator features
from ultra-conserved TF binding sites, and Task 2 / 3 by adding indicator features
from epigenetic regional data.

## Findings

First, we ran our dataset with the fasta files provided by Beer Labs and
obtained an average 5-fold cv score of about 0.82 for au-ROC.

We then ran on our own Vista Dataset and obtained around 0.85 with
count normalization for Task 1. For Task 2, our results ranged from 0.5 to
0.79. For task 3, the average au-ROC was 0.55.

Adding ultra-conserved indicator features boosted Task 1 to 0.88.
Adding epigenetic features boosted Task 2's average to around 0.82, and 0.77
for Task 3.

