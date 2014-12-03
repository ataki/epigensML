POS_INT="data/bed/pos_intersect.bed"
NEG_INT="data/bed/neg_intersect.bed"
HEART_INT="data/bed/heart_intersect.bed"
LIMB_INT="data/bed/limb_intersect.bed"
BRAIN_INT="data/bed/brain_intersect.bed"
FOREBRAIN_INT="data/bed/forebrain_intersect.bed"

POS="data/fasta/vt_hm_pos.fa"
NEG="data/fasta/vt_hm_neg.fa"
ANNOTATED="data/fasta/hm_annotated.fa"

TF="data/feature_bed/TF_binding_sites.bed"
HEART_MM="data/feature_bed/heart_mnemonics.bed"
BRAIN_MM="data/feature_bed/brain_mnemonics.bed"
LIMB_MM="data/feature_bed/limb_mnemonics.bed"

OUT_POS="out/pos_tf.npy"
OUT_NEG="out/neg_tf.npy"
OUT_HEART="out/heart.npy"
OUT_LIMB="out/limb.npy"
OUT_BRAIN="out/brain.npy"
OUT_FOREBRAIN="out/forebrain.npy"

TF_OP="tf"
NOT_TF_OP="not_tf"

#echo "POS TF"
#python src/gen_extra.py $POS_INT $POS $OUT_POS $TF_OP

#echo "NEG TF"
#python src/gen_extra.py $NEG_INT $NEG $OUT_NEG $TF_OP

#echo "HEART"
#python src/gen_extra.py $HEART_INT $ANNOTATED $OUT_HEART heart

#echo "LIMB"
#python src/gen_extra.py $LIMB_INT $ANNOTATED $OUT_LIMB limb

#echo "BRAIN"
#python src/gen_extra.py $BRAIN_INT $ANNOTATED $OUT_BRAIN brain

echo "Forebrain"
python src/gen_extra.py $FOREBRAIN_INT $ANNOTATED $OUT_FOREBRAIN forebrain
