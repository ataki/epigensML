#!/bin/bash

POS="data/fasta/hm_annotated.fa"
HEART="out/heart.npy"
LIMB="out/limb.npy"
BRAIN="out/brain.npy"

echo "Predicting tissues (limb|brain|heart|neural|other)"
echo "Using $POS"
echo "python src/tissue_clf.py $POS --heart=$HEART --limb=$LIMB --brain=$BRAIN"
python src/tissue_clf.py $POS --heart=$HEART --limb=$LIMB --brain=$BRAIN
