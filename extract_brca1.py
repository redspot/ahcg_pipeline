from __future__ import print_function
import csv
import io
from collections import defaultdict
from subprocess import check_output, CalledProcessError

coord_fn = '/home/basespace/hg19_refGene.txt'
reference_fn = '/home/basespace/resources/genome/hg19.fa'
bedtools_cmd = [
    'bedtools',
    'getfasta',
    '-fi',
    reference_fn,
    '-s',
    '-bed',
]
brca1_dict = defaultdict(dict)
DEBUG = True


def log_mesg(mesg):
    if DEBUG is True:
        print(mesg)

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
        if gene_name == 'BRCA1':
            log_mesg("found {} in {}".format(gene_id, coord_fn))
            brca1_dict[gene_id]['chrom'] = chrom
            brca1_dict[gene_id]['strand'] = strand
            brca1_dict[gene_id]['start'] = start
            brca1_dict[gene_id]['stop'] = stop
            brca1_dict[gene_id]['count'] = exon_count
            brca1_dict[gene_id]['starts'] = starts
            brca1_dict[gene_id]['stops'] = stops

            brca1_dict[gene_id]['prefix'] = start - starts[0]
            log_mesg("{}: start={} starts[0]={} prefix={}".format(
                gene_id,
                start,
                starts[0],
                brca1_dict[gene_id]['prefix']
            ))
            previous = 0
            for s, e in reversed(list(zip(starts, stops))):
                log_mesg("{}: s={} stop={}".format(
                    gene_id, s, stop,
                ))
                if stop > s:
                    brca1_dict[gene_id]['suffix'] = e - stop + previous
                    log_mesg("{}: e={} stop={} suffix={}".format(
                        gene_id, e, stop, brca1_dict[gene_id]['suffix']
                    ))
                    break
                else:
                    previous += e - s

for gene_id in brca1_dict.keys():
    bed_fn = gene_id + '.bed'
    starts = brca1_dict[gene_id]['starts']
    stops = brca1_dict[gene_id]['stops']
    #chrom,s,e,name,1,strand
    bed_fmt = "{}\t{}\t{}\t{}\t1\t{}\n"
    log_mesg("writing {}".format(bed_fn))
    with io.open(bed_fn, 'w') as bed_fd:
        for i, t in enumerate(zip(starts, stops)):
            s, e = t
            name = gene_id + '_' + str(i)
            bed_line = bed_fmt.format(
                brca1_dict[gene_id]['chrom'],
                s,
                e,
                name,
                brca1_dict[gene_id]['strand'],
            )
            bed_fd.write(bed_line)
    log_mesg("done writing {}".format(bed_fn))

    try:
        cmd = list(bedtools_cmd + [bed_fn])
        log_mesg("running {}".format(cmd))
        exons_fa = check_output(cmd, universal_newlines=True)
        log_mesg("done running: len={}".format(len(exons_fa)))
    except CalledProcessError:
        exons_fa = None
    if exons_fa is None:
        log_mesg("error: bedtools returned non-zero")
    else:
        fa_fn = gene_id + '.fa'
        with io.open(fa_fn, 'w') as fa_fd:
            fa_fd.write(exons_fa)
        exons = exons_fa.split('\n')[:-1]
        log_mesg("exons expected={}, exons got={}".format(
            brca1_dict[gene_id]['count'],
            len(exons) / 2
        ))
        sequence = ""
        prefix = brca1_dict[gene_id]['prefix']
        suffix = brca1_dict[gene_id]['suffix']
        if brca1_dict[gene_id]['strand'] == '-':
            exons = reversed(exons)
            new_suffix = prefix
            prefix = suffix
            suffix = new_suffix
        for line in exons:
            if not line.startswith('>'):
                sequence += line
        log_mesg("prefix={} suffix={} slice={}:{}".format(
            prefix,
            suffix,
            prefix,
            -suffix,
        ))
        sequence_trimmed = sequence[prefix:-suffix]
        log_mesg("first codon = {}".format(sequence_trimmed[0:3]))
        seq_fn = gene_id + '.seq.txt'
        log_mesg("writing sequence to {}".format(seq_fn))
        with io.open(seq_fn, 'w') as seq_fd:
            seq_fd.write(sequence_trimmed)
            seq_fd.write('\n')
            seq_fd.write(sequence)
