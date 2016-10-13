import sys
import csv
from collections import defaultdict
import vcf

cache = defaultdict(lambda: defaultdict(list))
patho = {}

for row in csv.DictReader(open(sys.argv[1])):
    patho_key = row['Genomic_Coordinate_hg37']
    if patho_key[0:3] == 'chr':
        patho_key = patho_key[3:]

    patho[patho_key] = row['Pathogenicity_research']
    recs = patho_key.split(':')
    cache[recs[0]][recs[1]].append(recs[2])

for rec in vcf.Reader(open(sys.argv[2])):
    for alt in rec.ALT:
        change = rec.REF + '>' + alt
        if change in cache[rec.CHROM][rec.POS]:
            path_key = ':'.join([rec.CHROM, rec.POS, change])
            if path_key in patho:
                print(rec, patho[path_key])
            else:
                print(rec, 'no pathogenicity info')
