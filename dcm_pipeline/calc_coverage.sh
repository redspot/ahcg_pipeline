#!/bin/bash
if [ "$#" -ne 2 ]
then
    echo "usage: $0 assembled.bam outdir"
    exit 1
fi
set -e

BASE=/home/basespace
GENE_BED="$BASE/resources/dcm/dcm_gene_list.bed"
GENELIST=$(cut -f4 $GENE_BED | sort -u | xargs)
IN_BAM="$1"
OUTDIR="$2"

samtools view -L "$GENE_BED" "$IN_BAM" -b > "$OUTDIR/dcm.bam"
bedtools genomecov -ibam "$OUTDIR/dcm.bam" -bga > "$OUTDIR/dcm_all.bga"
bedtools intersect -loj -split -a "$GENE_BED" -b "$OUTDIR/dcm_all.bga" -bed > "$OUTDIR/dcm_all.bed"
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$10,$6)}' "$OUTDIR/dcm_all.bed" > $OUTDIR/depths.bed
python $BASE/ahcg_pipeline/cov.py $OUTDIR/depths.bed $OUTDIR/depths.txt
test -s $OUTDIR/depths.txt && Rscript $BASE/ahcg_pipeline/draw_depth.R $OUTDIR/depths.txt $OUTDIR/depths.png || true
test -s $OUTDIR/depths.txt && Rscript $BASE/ahcg_pipeline/freq_plot.R $OUTDIR/depths.txt $OUTDIR/freq.png || true

for gene in $GENELIST
do
    awk '$4~/'${gene}'/ {print}' $GENE_BED > "$OUTDIR/${gene}.bed"
    samtools view -L "$OUTDIR/${gene}.bed" "$IN_BAM" -b > "$OUTDIR/${gene}.bam"
    bedtools genomecov -ibam "$OUTDIR/${gene}.bam" -bga > "$OUTDIR/${gene}.bga"
    bedtools intersect -loj -split -a "$OUTDIR/${gene}.bed" -b "$OUTDIR/${gene}.bga" -bed > "$OUTDIR/${gene}.bed"
    awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$10,$6)}' "$OUTDIR/${gene}.bed" > $OUTDIR/${gene}-depths.bed
    python $BASE/ahcg_pipeline/cov.py $OUTDIR/${gene}-depths.bed $OUTDIR/${gene}-depths.txt
    test -s $OUTDIR/${gene}-depths.txt && Rscript $BASE/ahcg_pipeline/draw_depth.R $OUTDIR/${gene}-depths.txt $OUTDIR/${gene}-depths.png || true
    test -s $OUTDIR/${gene}-depths.txt && Rscript $BASE/ahcg_pipeline/freq_plot.R $OUTDIR/${gene}-depths.txt $OUTDIR/${gene}-freq.png || true
done
