#!/bin/bash

POS="data/fasta/vt_hm_combo.fa"
NEG="data/fasta/vt_hm_neg.fa"

echo "Predicting combo regions vs non-combo regions"
echo "pos: $POS, neg: $NEG"
python src/enhancer_clf.py $POS $NEG

