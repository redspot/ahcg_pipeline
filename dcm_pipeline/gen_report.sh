#!/bin/bash
if [ "$#" -ne 2 ]
then
    echo "usage: $0 prefix outdir"
    exit 1
fi
set -e

BASE_DIR="/home/basespace"
GENE_BED="$BASE_DIR/resources/dcm/dcm_gene_list.bed"
GENELIST=$(cut -f4 $GENE_BED | sort -u | xargs)
TEMPLATES=$BASE_DIR/dcm_templates
PREFIX="$1"
OUTDIR="$2"
REPORT="$OUTDIR/${PREFIX}-final-report.txt"

test -s $OUTDIR/xref_pathogenic.csv || echo No Matches Found > $OUTDIR/xref-pathogenic.csv
test -s $OUTDIR/xref_benign.csv || echo No Matches Found > $OUTDIR/xref-benign.csv
test -s $OUTDIR/xref_conflicting.csv || echo No Matches Found > $OUTDIR/xref-conflicting.csv
test -s $OUTDIR/xref_other.csv || echo No Matches Found > $OUTDIR/xref-other.csv

cat \
<(echo -n $PREFIX) $TEMPLATES/article-title.txt \
$TEMPLATES/pathogenic-header.txt \
$TEMPLATES/table-header.txt \
$OUTDIR/xref_pathogenic.csv \
$TEMPLATES/table-footer.txt \
$TEMPLATES/benign-header.txt \
$TEMPLATES/table-header.txt \
$OUTDIR/xref_benign.csv \
$TEMPLATES/table-footer.txt \
$TEMPLATES/conflicting-header.txt \
$TEMPLATES/table-header.txt \
$OUTDIR/xref_conflicting.csv \
$TEMPLATES/table-footer.txt \
$TEMPLATES/other-header.txt \
$TEMPLATES/table-header.txt \
$OUTDIR/xref_other.csv \
$TEMPLATES/table-footer.txt \
<(sed dcm_templates/depth-chart.txt -e 's#XXXX#'$OUTDIR/depths.png'#') \
> $REPORT

for gene in $GENELIST
do
    if [ -s $OUTDIR/${gene}-depths.txt ]
    then
        cat <(sed dcm_templates/depth-chart.txt -e 's#XXXX#'$OUTDIR/${gene}-depths.png'#') \
        >> $REPORT
    fi
done

a2x -f pdf -d article \
--dblatex-opts '-P doc.layout="mainmatter" -P doc.publisher.show=0' \
$REPORT
