from __future__ import print_function
import csv
import sys
import io
from collections import defaultdict
#from subprocess import check_output, CalledProcessError

coord_fn = '/Users/wmartin45/school/biol8803f/hg19_refGene.txt'
reference_fn = '/Users/wmartin45/school/biol8803f/resources/genome/hg19.fa'
if len(sys.argv) < 2:
    gene_nums_fn = '/Users/wmartin45/school/biol8803f/bc_oc_NMnums.txt'
else:
    gene_nums_fn = sys.argv[1]

bedtools_cmd = [
    'bedtools',
    'getfasta',
    '-fi',
    reference_fn,
    '-s',
    '-fo',
    '-',
    '-bed',
]
bcoc_dict = defaultdict(dict)
gene_nums = []
DEBUG = True


def log_mesg(mesg):
    if DEBUG is True:
        print(mesg)

with io.open(gene_nums_fn) as gene_fd:
    gene_nums = [l.strip() for l in gene_fd.readlines()]

with io.open(coord_fn) as coord_fd:
    #1 = gene_id
    #2 = chrom
    #3 = strand
    #6 = coding start
    #7 = coding stop
    #8 = exon count
    #9 = exon starts
    #10 = exon stops
    #12 = gene name
    for row in csv.reader(coord_fd, delimiter='\t'):
        gene_id = row[1]
        chrom = row[2]
        strand = row[3]
        start = int(row[6])
        stop = int(row[7])
        exon_count = int(row[8])
        starts = list(map(int, filter(lambda s: s != "", row[9].split(','))))
        stops = list(map(int, filter(lambda s: s != "", row[10].split(','))))
        gene_name = row[12]
        if gene_id in gene_nums:
            bcoc_dict[gene_id]['chrom'] = chrom
            bcoc_dict[gene_id]['strand'] = strand
            bcoc_dict[gene_id]['start'] = start
            bcoc_dict[gene_id]['stop'] = stop
            bcoc_dict[gene_id]['count'] = exon_count
            bcoc_dict[gene_id]['starts'] = starts
            bcoc_dict[gene_id]['stops'] = stops
            bcoc_dict[gene_id]['name'] = gene_name

            bcoc_dict[gene_id]['prefix'] = start - starts[0]
            previous = 0
            for s, e in reversed(list(zip(starts, stops))):
                if stop > s:
                    bcoc_dict[gene_id]['suffix'] = e - stop + previous
                    break
                else:
                    previous += e - s

for gene_id in bcoc_dict.keys():
    bed_fn = gene_id + '.bed'
    start = bcoc_dict[gene_id]['start']
    stop = bcoc_dict[gene_id]['stop']
    #chrom,s,e,name,1,strand
    bed_fmt = "{}\t{}\t{}\t{}\t1\t{}\n"
    with io.open(bed_fn, 'w') as bed_fd:
        name = gene_id + ':' + bcoc_dict[gene_id]['name']
        bed_line = bed_fmt.format(
            bcoc_dict[gene_id]['chrom'],
            start,
            stop,
            name,
            bcoc_dict[gene_id]['strand'],
        )
        bed_fd.write(bed_line)
