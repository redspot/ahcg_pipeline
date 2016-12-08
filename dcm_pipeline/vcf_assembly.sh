#!/bin/bash
if [ "$#" -ne 3 ]
then
    echo "usage: $0 reads1.fastq reads2.fastq outdir"
    exit 1
fi
set -e

BASE_DIR="/home/basespace"
TRIM=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTER=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
BOWTIE=/home/basespace/.linuxbrew/bin/bowtie2
PICARD=$BASE_DIR/ahcg_pipeline/lib/picard.jar
GATK=$BASE_DIR/ahcg_pipeline/lib/GenomeAnalysisTK.jar

python3 $BASE_DIR/ahcg_pipeline/ahcg_pipeline.py \
    -t $TRIM \
    -b $BOWTIE \
    -p $PICARD \
    -g $GATK \
    -i $1 $2 \
    -w $BASE_DIR/resources/genome/hg19 \
    -d $BASE_DIR/resources/dbsnp/dbsnp_138.hg19.vcf.gz \
    -r $BASE_DIR/resources/genome/hg19.fa \
    -a $ADAPTER \
    -o $3

#samtools faidx hg19.fa
#picard CreateSequenceDictionary R=hg19.fa O=hg19.dict
#bowtie2-build hg19.fa hg19
#tabix -p vcf dbsnp_138.hg19.vcf.gz
