#!/bin/bash
if [ "$#" -eq 0 ]
then
    echo "usage: $0 variants_from_bam.vcf"
    exit 1
fi

BASE=/Users/wmartin45/school/biol8803f
JAR=$BASE/ahcg_pipeline/lib/GenomeAnalysisTK.jar
REF=$BASE/resources/genome/hg19.fa
#VARIANTS=$BASE/hw6/patient2_variants.vcf
VARIANTS="$1"
RECAL=$(dirname $VARIANTS)/hw6/output.recal
TRANCH=$(dirname $VARIANTS)/hw6/output.tranches
java -jar $JAR \
    -T ApplyRecalibration \
    -R $REF \
    -input $VARIANTS \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile $RECAL \
    --tranches_file $TRANCH \
    -o NA12878_variants_recal.vcf
