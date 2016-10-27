from __future__ import print_function
import vcf
from io import StringIO
import sys

if len(sys.argv) < 3:
    print("usage: {} input.txt output.vcf".format(sys.argv[0]))
    exit(1)

header_txt = u"""##fileformat=VCFv4.1
##INFO=<ID=PG,Number=1,Type=String,Description="Pathogenicity">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
"""

header_fd = StringIO(header_txt)
reader = vcf.Reader(header_fd)

with open(sys.argv[1]) as in_fd,\
        open(sys.argv[2], 'w') as out_fd:
    writer = None
    try:
        writer = vcf.Writer(out_fd, reader)
        for line in in_fd:
            change, desc = line.split('\t')
            desc = desc.strip().replace(";", "|")
            chrom, pos, snp = change.split(':')
            pos = int(pos)
            ref, alt = snp.split('>')
            alt = [vcf.model._Substitution(alt)]
            info = {"PG": desc}
            rec = vcf.model._Record(chrom, pos, None, ref, alt, 1.0, [], info, ".", {})
            writer.write_record(rec)
    finally:
        if writer:
            writer.close()

#CHROM str
#POS int
#ID None
#REF str
#ALT list of vcf.model._Substitution
#QUAL float
#FILTER empty list
#INFO dict
#FORMAT str
#sample_indexes dict
