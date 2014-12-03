#!/bin/bash

POS="data/fasta/hm_annotated.fa"

echo "Predicting tissues (fore|mid|hind)"
echo "Using $POS"
python src/fg_clf.py $POS --prediction_type=binary

