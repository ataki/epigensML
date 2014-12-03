#!/bin/bash

POS="data/fasta/vt_hm_pos.fa"
NEG="data/fasta/vt_hm_neg.fa"

echo "Predicting general enhancer activity"
echo "pos: $POS, neg: $NEG"
python src/enhancer_clf.py $POS $NEG

