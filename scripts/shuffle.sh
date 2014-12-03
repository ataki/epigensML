#!/bin/bash

echo "Shuffling human genes"
python src/gen_neg_seqs.py data/bed/vt_hm.bed data/bed/vt_hm_shf.bed

echo "Shuffling mouse genes"
python src/gen_neg_seqs.py data/bed/vt_mm.bed data/bed/vt_mm_shf.bed mm9

