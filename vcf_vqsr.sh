#!/bin/bash
if [ "$#" -ne 2 ]
then
    echo "usage: $0 assembled.vcf outdir"
    exit 1
fi
set -e

BASE=/home/basespace
JAR=$BASE/ahcg_pipeline/lib/GenomeAnalysisTK.jar
REF=$BASE/resources/genome/hg19.fa
VARIANTS="$1"
OUTBASE="$2"
HAPMAP=$BASE/resources/vqsr/hapmap_3.3.hg19.sites.vcf.gz
OMNI=$BASE/resources/vqsr/1000G_omni2.5.hg19.sites.vcf.gz
PHASE=$BASE/resources/vqsr/1000G_phase1.snps.high_conf.hg19.sites.vcf.gz
DBSNP=$BASE/resources/dbsnp/dbsnp_138.hg19.vcf.gz
java -Xmx4g -jar $JAR \
	-T VariantRecalibrator \
	-R $REF \
	-input $VARIANTS \
	-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
	-resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI \
	-resource:1000G,known=false,training=true,truth=false,prior=10.0 $PHASE \
	-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP \
	-an DP \
	-an QD \
	-an FS \
	-an SOR \
	-an MQ \
	-an MQRankSum \
	-an ReadPosRankSum \
	-mode SNP \
	-recalFile "$OUTBASE/output.recal" \
	-tranchesFile "$OUTBASE/output.tranches"
