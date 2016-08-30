#!/bin/bash

#bail out if there are any errors
set -e

#setup variables
TRIM_PATH=
BOWTIE_PATH=
PICARD_PATH=
GATK_PATH=
declare -a INPUT_FILES
REFERENCE=
BT2_INDICES=
DBSNP=
ADAPTER=
OUTPUT_PATH=
TRIM_SUFFIX=".trimmed"
UNUSED_SUFFIX=".unused"
SAM_PATH=
READGROUP_PATH=
MARKDUP_PATH=
METRICS=
FIXMATE_PATH=
INTERVALS=
INDEL_PATH=
TABLE_PATH=
FINAL_PATH=
VARIANTS=

#parse args
while getopts "t:b:p:g:w:i:d:r:a:o:" opt; do
    case $opt in
        t) TRIM_PATH="$OPTARG" ;;
        b) BOWTIE_PATH="$OPTARG" ;;
        p) PICARD_PATH="$OPTARG" ;;
        g) GATK_PATH="$OPTARG" ;;
        w) BT2_INDICES="$OPTARG" ;;
        i) INPUT_FILES+=("$OPTARG") ;;
        d) DBSNP="$OPTARG" ;;
        r) REFERENCE="$OPTARG" ;;
        a) ADAPTER="$OPTARG" ;;
        o) OUTPUT_PATH="$OPTARG" ;;
    esac
done

#check paths
[ -f "$" ] || (echo "file not found: $"; exit 1)
[ -f "$TRIM_PATH" ] || (echo "file not found: $TRIM_PATH"; exit 1)
[ -f "$BOWTIE_PATH" ] || (echo "file not found: $BOWTIE_PATH"; exit 1)
[ -f "$PICARD_PATH" ] || (echo "file not found: $PICARD_PATH"; exit 1)
[ -f "$GATK_PATH" ] || (echo "file not found: $GATK_PATH"; exit 1)
[ -z "$(echo ${BT2_INDICES}.*.bt2)" ] && (echo "bowtie2 indices not found: $BT2_INDICES"; exit 1)
for val in "${INPUT_FILES[@]}"; do
    [ -f "$val" ] || (echo "file not found: $val"; exit 1)
done
[ "${#INPUT_FILES[@]}" -eq "2" ] || (echo "expecting exactly 2 input files"; exit 1)
[ -f "$DBSNP" ] || (echo "file not found: $DBSNP"; exit 1)
[ -f "$REFERENCE" ] || (echo "file not found: $REFERENCE"; exit 1)
[ -f "$ADAPTER" ] || (echo "file not found: $ADAPTER"; exit 1)

#create output path
[ -d "$OUTPUT_PATH" ] || mkdir -p "$OUTPUT_PATH"

#define additional paths
BASENAME1="$(basename ${INPUT_FILES[0]})"
BASENAME2="$(basename ${INPUT_FILES[1]})"
TRIM_NAME1="${BASENAME1}${TRIM_SUFFIX}"
TRIM_NAME2="${BASENAME2}${TRIM_SUFFIX}"
UNUSED_NAME1="${BASENAME1}${UNUSED_SUFFIX}"
UNUSED_NAME2="${BASENAME2}${UNUSED_SUFFIX}"
SAM_PATH="${OUTPUT_PATH}/${TRIM_NAME1}.sam"
READGROUP_PATH="${SAM_PATH}.readgroup.bam"
MARKDUP_PATH="${SAM_PATH}.markdup.bam"
METRICS="${SAM_PATH}.markdup.metrics"
FIXMATE_PATH="${SAM_PATH}.fixmate.bam"
INTERVALS="${SAM_PATH}.intervals"
INDEL_PATH="${SAM_PATH}.indel.bam"
TABLE_PATH="${SAM_PATH}.table"
FINAL_PATH="${SAM_PATH}.final.bam"
VARIANTS="${OUTPUT_PATH}/variants.vcf"

#run trim
java -jar "$TRIM_PATH" \
 PE \
 -phred33 \
 "${INPUT_FILES[0]}" \
 "${INPUT_FILES[1]}" \
 "${INPUT_FILES[0]}$TRIM_SUFFIX" \
 "${INPUT_FILES[0]}$UNUSED_SUFFIX" \
 "${INPUT_FILES[1]}$TRIM_SUFFIX" \
 "${INPUT_FILES[1]}$UNUSED_SUFFIX" \
 "ILLUMINACLIP:$ADAPTER:2:30:10" \
 LEADING:0 \
 TRAILING:0 \
 SLIDINGWINDOW:4:15 \
 MINLEN:36

#align using bowtie2
/usr/local/bin/bowtie2 \
 -x "$BT2_INDICES" \
 -S "$SAM_PATH" \
 -p 1 \
 -1 "${INPUT_FILES[0]}$TRIM_SUFFIX" \
 -2 "${INPUT_FILES[1]}$TRIM_SUFFIX"

#Add read group information
java -Xmx1g -jar "$PICARD_PATH" \
 AddOrReplaceReadGroups \
 "I=$SAM_PATH" \
 "O=${READGROUP_PATH}" \
 SORT_ORDER=coordinate \
 RGID=Test \
 RGLB=ExomeSeq \
 RGPL=Illumina \
 RGPU=HiSeq2500 \
 RGSM=Test \
 RGCN=AtlantaGenomeCenter \
 RGDS=ExomeSeq \
 RGDT=2016-08-24 \
 RGPI=null \
 RGPG=Test \
 RGPM=Test \
 CREATE_INDEX=true

#Mark PCR duplicates
java -Xmx1g -jar "$PICARD_PATH" \
 MarkDuplicates \
 "I=${READGROUP_PATH}" \
 "O=${MARKDUP_PATH}" \
 "METRICS_FILE=${METRICS}" \
 REMOVE_DUPLICATES=false \
 ASSUME_SORTED=true \
 CREATE_INDEX=true

#Fix mate information
java -Xmx1g -jar "$PICARD_PATH" \
 FixMateInformation \
 "I=${MARKDUP_PATH}" \
 "O=${FIXMATE_PATH}" \
 ASSUME_SORTED=true \
 ADD_MATE_CIGAR=true \
 CREATE_INDEX=true

#Run realigner target creator
java -jar "$GATK_PATH" \
 -T RealignerTargetCreator \
 -o "$INTERVALS" \
 -nt 1 \
 -I "${FIXMATE_PATH}" \
 -R "$REFERENCE" \
 -known "$DBSNP"

#Run indel realigner
java -jar "$GATK_PATH" \
 -T IndelRealigner \
 --targetIntervals "$INTERVALS" \
 -o "${INDEL_PATH}" \
 -I "${FIXMATE_PATH}" \
 -R "$REFERENCE"

#Base quality score recalibration
java -jar "$GATK_PATH" \
 -T BaseRecalibrator \
 -R "$REFERENCE" \
 -I "${INDEL_PATH}" \
 -o "${TABLE_PATH}" \
 -nct 1 \
 -cov ReadGroupCovariate \
 -knownSites "$DBSNP"

#Print Reads
java -jar "$GATK_PATH" \
 -T PrintReads \
 -R "$REFERENCE" \
 -I "${INDEL_PATH}" \
 -o "${FINAL_PATH}" \
 -BQSR "${TABLE_PATH}" \
 -nct 1

#Haplotype caller
java -jar "$GATK_PATH"
 -T HaplotypeCaller \
 -R "$REFERENCE" \
 -I "${FINAL_PATH}" \
 --dbsnp "$DBSNP" \
 -o "$VARIANTS"
 -nct 1 \
 -gt_mode DISCOVERY

echo Variant call pipeline completed
echo VCF file can be found at $VARIANTS
exit
test_r1.fastq_trimmed.fq
test_r1.fastq_trimmed.sam
test_r1.fastq_unused.fq
test_r2.fastq_trimmed.fq
test_r2.fastq_unused.fq
