#!/bin/bash
if [ "$#" -ne 4 ]
then
    echo "usage: $0 in_recal.vcf out_dcm.vcf out_xref.vcf outdir"
    exit 1
fi
set -e

BASE=/home/basespace
GENE_BED="$BASE/resources/dcm/dcm_gene_list.bed"
CLN_DCM="$BASE/resources/dcm/clinvar_dcm_genes.vcf"
IN_VCF="$1"
DCM_VCF="$2"
XREF_VCF="$3"
OUTDIR="$4"

#strip chromosome from variants lines
sed -i.bak -E -e 's/^chr//' "$IN_VCF"

#shrink variants to just DCM genes
bedtools intersect -a "$IN_VCF" -b "$GENE_BED" -header > "$DCM_VCF"

#match variants to clinvar
bedtools intersect -b "$DCM_VCF" -a "$CLN_DCM" -header > "$XREF_VCF"

#generate simple report on findings
python3 $BASE/ahcg_pipeline/parse_clnsig.py -i "$XREF_VCF" \
-p "$OUTDIR/xref_pathogenic.txt" \
-b "$OUTDIR/xref_benign.txt" \
-c "$OUTDIR/xref_conflicting.txt" \
-o "$OUTDIR/xref_other.txt"

csvtool -t TAB -u COMMA cat "$OUTDIR/xref_pathogenic.txt" > "$OUTDIR/xref_pathogenic.csv"
csvtool -t TAB -u COMMA cat "$OUTDIR/xref_benign.txt" > "$OUTDIR/xref_benign.csv"
csvtool -t TAB -u COMMA cat "$OUTDIR/xref_conflicting.txt" > "$OUTDIR/xref_conflicting.csv"
csvtool -t TAB -u COMMA cat "$OUTDIR/xref_other.txt" > "$OUTDIR/xref_other.csv"
