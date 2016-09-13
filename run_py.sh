#!/bin/bash
BASE_DIR="/Users/wmartin45/school/biol8803f"
TRIM=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTER=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
BOWTIE=/usr/local/bin/bowtie2
PICARD=$BASE_DIR/ahcg_pipeline/lib/picard.jar
GATK=$BASE_DIR/ahcg_pipeline/lib/GenomeAnalysisTK.jar

python3 $BASE_DIR/ahcg_pipeline/ahcg_pipeline.py \
    -t $TRIM \
    -b $BOWTIE \
    -p $PICARD \
    -g $GATK \
    -i $BASE_DIR/reads_B1_27000x150bp_0S_0I_0D_0U_0N_1.fq.gz $BASE_DIR/reads_B1_27000x150bp_0S_0I_0D_0U_0N_2.fq.gz \
    -w $BASE_DIR/resources/genome/chr17 \
    -d $BASE_DIR/resources/dbsnp/dbsnp_138.hg19.vcf.gz \
    -r $BASE_DIR/resources/genome/chr17.fa \
    -a $ADAPTER \
    -o $BASE_DIR/hw1

#samtools faidx hg19.fa
#picard CreateSequenceDictionary R=hg19.fa O=hg19.dict
#bowtie2-build hg19.fa hg19
#tabix -p vcf dbsnp_138.hg19.vcf.gz
