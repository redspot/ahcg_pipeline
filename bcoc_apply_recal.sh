#!/bin/bash
BASE=/Users/wmartin45/school/biol8803f
JAR=$BASE/ahcg_pipeline/lib/GenomeAnalysisTK.jar
REF=$BASE/resources/genome/hg19.fa
VARIANTS=$BASE/hw3/NA12878_variants.vcf
RECAL=$BASE/hw3/output.recal
TRANCH=$BASE/hw3/output.tranches
java -jar $JAR \
    -T ApplyRecalibration \
    -R $REF \
    -input $VARIANTS \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile $RECAL \
    --tranches_file $TRANCH \
    -o NA12878_variants_recal.vcf
