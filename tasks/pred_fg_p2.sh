#!/bin/bash

POS="data/fasta/hm_annotated.fa"
FOREBRAIN="out/forebrain.npy"

echo "Predicting tissues (fore|mid|hind)"
echo "Using $POS"
python src/fg_clf.py $POS --forebrain=$FOREBRAIN --prediction_type=binary

