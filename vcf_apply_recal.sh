#!/bin/bash
if [ "$#" -ne 3 ]
then
    echo "usage: $0 assembled.vcf recal_output.vcf recaldir"
    exit 1
fi
set -e

BASE=/home/basespace
JAR=$BASE/ahcg_pipeline/lib/GenomeAnalysisTK.jar
REF=$BASE/resources/genome/hg19.fa
VARIANTS="$1"
RECAL_OUT="$2"
OUTBASE="$3"
RECAL="$OUTBASE/output.recal"
TRANCH="$OUTBASE/output.tranches"
java -jar $JAR \
    -T ApplyRecalibration \
    -R $REF \
    -input $VARIANTS \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile $RECAL \
    --tranches_file $TRANCH \
    -o $RECAL_OUT
