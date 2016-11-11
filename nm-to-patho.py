#!/usr/bin/env python3
"""links to Entrez API
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc110
https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Understanding_the_Eutilities_Wi
https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ELink
https://eutils.ncbi.nlm.nih.gov/corehtml/query/static/entrezlinks.html
https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
"""
import logging
#import json
import gzip
from io import StringIO
from xml.parsers.expat import ExpatError
#from xml.dom import minidom

import click

import vcf
from Bio import Entrez
from Bio.Entrez.Parser import (NotXMLError, CorruptedXMLError, ValidationError)

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)

Entrez.email = "wilson.martin@gtri.gatech.edu"
RETMAX = 100000  # Entrez allows a max of 100k returned


def check_handle(handle):
    data = "!!NotSet!!"
    try:
        #record = Entrez.read(handle)
        data = handle.read()
        record = Entrez.read(StringIO(data))
    except (NotXMLError, CorruptedXMLError, ValidationError, ExpatError):
        logger.exception("the results are corrupted:")
        logger.debug(data)
        exit(1)
    #except NotImplementedError:
    #    record = minidom.parseString(data)
    except RuntimeError as e:
        if e.args and 'Empty' in e.args[0]:
            logger.error("no results were returned. bailing out.")
            exit(1)
        else:
            raise
    except:
        logger.exception("something broke:")
        logger.debug(data)
        exit(1)
    total = 0
    if isinstance(record, Entrez.Parser.ListElement):
        for rec in record:
            total += 1 if isinstance(rec, str) else len(rec['IdList'])
    else:
        total = int(record['Count'])
    if total == 0:
        logger.error("no results were returned. bailing out.")
        exit(1)
    return record


@click.command()
@click.option('-l', '--locus', 'locus_fd', required=True, type=click.File('r'), help='path to tsv file of locus, 2nd column')
#@click.option('-j', '--json', 'json_fd', type=click.File('w'), help='output filename for fetched records as json')
@click.option('-c', '--clinvar', 'clinvar_fd', required=True, type=click.File('r'), help='path to clinvar vcf')
@click.option('-f', '--filtered', 'filt_fd', required=True, type=click.File('w'), help='output filename for filtered vcf')
@click.option('-z', '--compress', 'compress', is_flag=True, help='gzip compress output vcf?')
def main(
    locus_fd,
    #json_fd,
    clinvar_fd,
    filt_fd,
    compress,
):
    locus_lines = locus_fd.readlines()
    locus_list = [l.split('\t')[1].strip() for l in locus_lines]

    nuc_query = ' OR '.join([l + '[Accession]' for l in locus_list])

    history = []
    logger.info("doing esearch for {}".format(nuc_query))
    handle = Entrez.esearch(
        db="nuccore",
        term=nuc_query,
        usehistory="y",
    )
    record = check_handle(handle)
    history.append([
        record["WebEnv"],
        record["QueryKey"]
    ])
    logger.info("got {} entries".format(record['Count']))
    logger.info("got history {}".format(history[0]))

    logger.info("doing elink from nuccore to snp")
    handle = Entrez.elink(
        db="snp", dbfrom="nuccore",
        webenv=history[0][0],
        query_key=history[0][1],
        cmd="neighbor_history",
        linkname="nuccore_snp",
    )
    record = check_handle(handle)
    record = record[0]
    history.append([
        record["WebEnv"],
        record['LinkSetDbHistory'][0]["QueryKey"]
    ])
    logger.info("got {} entries".format(len(record['IdList'])))
    logger.info("got history {}".format(history[1]))

    logger.info("doing esearch to snp and filtering by clinical significance")
    #clinsig_filter = '(("clinsig pathogenic"[Properties]))'
    clinsig_filter = '"pathogenic snp"[Filter]'
    #https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Understanding_the_Eutilities_Wi
    # set term="#1 AND stuff[THING]" where '1' is query_key
    # then pass webenv and term to esearch, but not query_key
    snp_term = '#{} AND {}'.format(history[1][1], clinsig_filter)
    logger.debug("snp_term = '{}'".format(snp_term))
    handle = Entrez.esearch(
        db="snp",
        retmax=RETMAX,  # Entrez allows a max of 100k returned
        term=snp_term,
        webenv=history[1][0],
        usehistory="y",
    )
    record = check_handle(handle)
    history.append([
        record["WebEnv"],
        record["QueryKey"]
    ])
    logger.info("got {} entries".format(record['Count']))
    logger.info("got history {}".format(history[2]))

    #logger.info("doing efetch from snp")
    #handle = Entrez.efetch(
    #    db="snp",
    #    retmax=RETMAX,  # Entrez allows a max of 100k returned
    #    #retmode='json',
    #    #retmode='text',
    #    retmode='xml',
    #    #rettype='uilist',
    #    webenv=history[2][0],
    #    query_key=history[2][1],
    #    usehistory="y",
    #)
    #record = check_handle(handle)
    #logger.info("got {} entries".format(len(record)))

    try:
        count = int(record['Count'])
        #count = len(record['IdList'])
        #count = len(record)
    except:
        logger.error("cannot determine how many records were returned")
        logger.error("cowardly refusing to continue")
        exit(1)

    if count > RETMAX:
        logger.warn("only retrieving the first {} entries of {}".format(
            RETMAX, count))

    #if json_fd is not None:
    #    efetch_json = json.dumps(record)
    #    logger.info("writing json file, {} bytes".format(len(efetch_json)))
    #    json_fd.write(efetch_json)

    snp_idlist = record["IdList"]
    logger.info("filtering clinvar vcf by {} SNP ID's".format(len(snp_idlist)))
    try:
        clinvar_reader = vcf.Reader(clinvar_fd)
    except:
        clinvar_fn = clinvar_fd.name
        clinvar_fd.close()
        vcf_str = gzip.open(clinvar_fn, 'rt').read()
        clinvar_reader = vcf.Reader(StringIO(vcf_str))
    if compress:
        filt_fn = filt_fd.name
        filt_fd.close()
        filt_fd = gzip.open(filt_fn, 'wt')
    filt_writer = vcf.Writer(filt_fd, clinvar_reader)
    f_count = 0
    for i, row in enumerate(clinvar_reader):
        rs = str(row.INFO['RS'])
        if rs in snp_idlist:
            f_count += 1
            filt_writer.write_record(row)
    if compress:
        filt_fd.close()  # explicit close() just in case it's gzip'd
    logger.info("filtered down to {} entries from {} total".format(f_count, i + 1))

if __name__ == '__main__':
    main()
