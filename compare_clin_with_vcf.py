import sys
import csv
from collections import defaultdict
import vcf

cache = defaultdict(lambda: defaultdict(list))
patho = {}

if len(sys.argv) < 3:
    print("usage: {} file.vcf file1.csv file2.csv ...".format(sys.argv[0]))
    exit(1)

for fn in sys.argv[2:]:
    for row in csv.DictReader(open(fn)):
        patho_key = row['Genomic_Coordinate_hg37']
        if patho_key[0:3] == 'chr':
            patho_key = patho_key[3:]
        if row['Pathogenicity_research'] != "":
            patho[patho_key] = row['Pathogenicity_research']
        recs = patho_key.split(':')
        cache[recs[0]][recs[1]].append(recs[2])

for rec in vcf.Reader(open(sys.argv[1])):
    chrom = rec.CHROM[3:] if rec.CHROM[0:3] == 'chr' else str(rec.CHROM)
    if chrom in cache:
        pos = str(rec.POS)
        for alt in map(str, rec.ALT):
            change = rec.REF + '>' + alt
            if change in cache[chrom][pos]:
                path_key = ':'.join([chrom, pos, change])
                if path_key in patho and patho[path_key]:
                    print(path_key + '\t' + patho[path_key])
