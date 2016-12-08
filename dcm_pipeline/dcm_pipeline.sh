#!/bin/bash
if [ "$#" -ne 4 ]
then
    echo "usage: $0 reads1.fastq reads2.fastq outdir prefix"
    exit 1
fi

set -e

BASE=/home/basespace

FQ1="$1"
FQ2="$2"
OUTDIR="$3"
PREFIX="$4"

fq1_base=$(basename "$FQ1")
fq1_dir=$(dirname "$FQ1")
fq1_ext="${fq1_base##*.}"
fq1_fn="${fq1_base%.*}"
aligned_bam="${OUTDIR}/${fq1_base}_trimmed_final.bam"
cleaned_bam="${OUTDIR}/${PREFIX}_cleaned.bam"
aligned_vcf="${OUTDIR}/variants.vcf"
recal_vcf="${OUTDIR}/${PREFIX}_recal.vcf"
dcm_vcf="${OUTDIR}/${PREFIX}_dcm.vcf"
xref_vcf="${OUTDIR}/${PREFIX}_xref.vcf"

./vcf_assembly.sh "$FQ1" "$FQ2" "$OUTDIR"

samtools view -H "$aligned_bam" \
| sed -e '/^@SQ/s/SN\:chr/SN\:/' \
| samtools reheader - "$aligned_bam" \
> "$cleaned_bam"

$BASE/ahcg_pipeline/vcf_vqsr.sh "$aligned_vcf" "$OUTDIR"
$BASE/ahcg_pipeline/vcf_apply_recal.sh "$aligned_vcf" "$recal_vcf" "$OUTDIR"
./xref_clinvar_dcm.sh "$recal_vcf" "$dcm_vcf" "$xref_vcf" "$OUTDIR" 
./calc_coverage.sh "$cleaned_bam" "$OUTDIR"
./gen_report "$PREFIX" "$OUTDIR"
