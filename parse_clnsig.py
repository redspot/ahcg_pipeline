#!/usr/bin/env python3
"""
parsing rules:
https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
"""
import logging
import gzip
from io import StringIO

import click

import vcf

logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s %(message)s',
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logger = logging.getLogger(__name__)


def get_clnset(c_list):
    return set(
        map(
            int,
            sum(
                [
                    cs.split('|')
                    for cs in c_list
                ],
                []
            )
        )
    )

#Pathogenic or Likely pathogenic
any_patho = set([4, 5])
#Benign or Likely benign
any_benign = set([2, 3])
#Pathogenic or Likely pathogenic or Benign or Likely benign
any_signif = any_patho | any_benign
#Pathogenic or Likely pathogenic or Benign or Likely benign or Uncertain
any_acmg = set([0]) | any_signif

descriptions = {
    0: 'Uncertain significance',
    1: 'Untested',
    2: 'Benign',
    3: 'Likely benign',
    4: 'Likely pathogenic',
    5: 'Pathogenic',
    6: 'Drug response',
    253: 'Conflicting data from submitters',
    254: 'Conflicting interpretations of pathogenicity',
    255: 'Risk factor',
}

categories = [
    'Pathogenic',
    'Benign',
    'Conflicting interpretations',
    'Non-pathogenic, non-benign, and non-conflicting',
    'Null, Impossible',
]


def sort_and_label_clnset(cs):
    patho_values = cs & any_patho
    benign_values = cs & any_benign
    signif_values = cs & any_signif
    acmg_values = cs & any_acmg
    non_acmg_values = cs - any_acmg
    is_uncertain = 0 in cs
    is_conflicting = False

    desc_list = []
    category = 4

    #(Pathogenic or Likely pathogenic or Benign or Likely benign) AND Uncertain significance
    if signif_values and is_uncertain:
        is_conflicting = True
    #(Pathogenic or Likely pathogenic) AND (Benign or Likely benign)
    if patho_values and benign_values:
        is_conflicting = True

    #Any single ACMG value AND any non-ACMG value
    if len(acmg_values) == 1 and len(non_acmg_values) == 1:
        val = acmg_values.pop()
        desc_list = [
            descriptions[val],
            descriptions[non_acmg_values.pop()],
        ]
        if val in any_patho:
            category = 0
        if val in any_benign:
            category = 1
    #Conflicting ACMG values AND any non-ACMG value
    elif is_conflicting and len(non_acmg_values) == 1:
        category = 2
        desc_list = [descriptions[254]]
        if is_uncertain:
            desc_list.append(descriptions[0])
        desc_list.append(descriptions[non_acmg_values.pop()])
    #No ACMG value AND multiple non-ACMG values
    elif len(acmg_values) == 0 and len(non_acmg_values) > 1:
        category = 3
        for v in non_acmg_values:
            desc_list.append(descriptions[v])
    #Any Conflicting ACMG values
    elif is_conflicting:
        category = 2
        desc_list = [descriptions[254]]
        if is_uncertain:
            desc_list.append(descriptions[0])
    #Any Benign
    elif benign_values:
        category = 1
        for v in benign_values:
            desc_list.append(descriptions[v])
    #Any Pathogenic
    elif patho_values:
        category = 0
        for v in sorted(patho_values, reverse=True):
            desc_list.append(descriptions[v])
    #Just Uncertain
    elif is_uncertain:
        category = 3
        desc_list = [descriptions[0]]
    #Something else
    else:
        category = 3
        desc_list = [descriptions[253]]

    description = '/'.join(desc_list)
    return (category, description)


@click.command()
@click.option('-i', '--input-vcf', 'in_fd', required=True, type=click.File('r'), help='path to input vcf with CLNSIG field')
def main(
    in_fd
):
    try:
        in_reader = vcf.Reader(in_fd)
    except:
        in_fn = in_fd.name
        in_fd.close()
        vcf_str = gzip.open(in_fn, 'rt').read()
        in_reader = vcf.Reader(StringIO(vcf_str))

    cln_sets = [[], [], [], []]

    for rec in in_reader:
        clnsig = rec.INFO['CLNSIG']
        rs = 'rs' + str(rec.INFO['RS'])
        gene = rec.INFO['GENEINFO']
        clnset = get_clnset(clnsig)
        cat, label = sort_and_label_clnset(clnset)
        cln_sets[cat].append((rs, gene, label))

    for i, cat in enumerate(cln_sets):
        logger.info('category: {}'.format(categories[i]))

        for rs, gene, label in cat:
            logger.info('|{:<15}|{:<15}|{}|'.format(rs, gene, label))

        logger.info('total: {}\n'.format(len(cat)))


if __name__ == '__main__':
    main()
