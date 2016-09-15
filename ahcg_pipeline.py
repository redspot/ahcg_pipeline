#!/usr/bin/env python3
from __future__ import print_function
import os
import sys
import glob
#import logging
import argparse
import subprocess
try:
    FileNotFoundError
except NameError:
    FileNotFoundError = IOError

THREAD_COUNT = str(8)


def run_process(args, mesg):
    p = subprocess.Popen(args, shell=False)
    p.wait()
    if p.returncode != 0:
        print('{0}; Exiting program'.format(mesg))
        sys.exit()


def main(trim_path, bowtie_path, picard_path, gatk_path,
         input_path, index_path, dbsnp_path, adapter_path,
         ref_path, out_path):

    #Get complete path
    trim_path = os.path.abspath(trim_path)
    bowtie_path = os.path.abspath(bowtie_path)
    picard_path = os.path.abspath(picard_path)
    gatk_path = os.path.abspath(gatk_path)
    input_path = [os.path.abspath(files) for files in input_path]
    index_path = os.path.abspath(index_path)
    dbsnp_path = os.path.abspath(dbsnp_path)
    ref_path = os.path.abspath(ref_path)
    out_path = os.path.abspath(out_path)
    adapter_path = os.path.abspath(adapter_path)

    #Check if paths exist
    path_checks = [
        (trim_path, 'Trimmomatic'),
        (bowtie_path, 'Bowtie'),
        (picard_path, 'Picard'),
        (gatk_path, 'Gatk'),
        (ref_path, 'Reference fasta file'),
        (dbsnp_path, 'dbSNP file'),
        (adapter_path, 'Adapter file'),
    ]
    for path, mesg in path_checks:
        if not os.path.exists(path):
            raise FileNotFoundError('{} not found at {}'.format(path, mesg))

    for files in input_path:
        if not os.path.exists(files):
            raise FileNotFoundError('Fastq files not found at {0}'.format(files))

    indicies = glob.glob('{0}.*.bt2'.format(index_path))
    #print(indicies)
    if len(indicies) == 0:
        raise FileNotFoundError('Bowtie index not found at {0}'.format(index_path))

    #Creat output directory
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    #build all of the commands
    read1 = input_path[0]
    read2 = input_path[1]
    tread1 = '{1}_trimmed.fq'.format(out_path, os.path.splitext(read1)[0])
    tread2 = '{1}_trimmed.fq'.format(out_path, os.path.splitext(read2)[0])
    sread1 = '{1}_unused.fq'.format(out_path, os.path.splitext(read1)[0])
    sread2 = '{1}_unused.fq'.format(out_path, os.path.splitext(read2)[0])
    tcmd = ['java', '-jar', trim_path, 'PE', '-phred33', read1, read2, tread1,
            sread1, tread2, sread2, 'ILLUMINACLIP:{0}:2:30:10'.format(adapter_path),
            'LEADING:0', 'TRAILING:0', 'SLIDINGWINDOW:4:15', 'MINLEN:36']
    sam_path = '{1}.sam'.format(out_path, os.path.splitext(tread1)[0])
    bcmd = [bowtie_path, '-x', index_path, '-S', sam_path, '-p', '1', '-1',
            tread1, '-2', tread2]
    add_path = '{0}/{1}_RG.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    acmd = ['java', '-Xmx1g', '-jar', picard_path, 'AddOrReplaceReadGroups',
            'I=' + sam_path, 'O=' + add_path, 'SORT_ORDER=coordinate', 'RGID=Test',
            'RGLB=ExomeSeq', 'RGPL=Illumina', 'RGPU=HiSeq2500', 'RGSM=Test',
            'RGCN=AtlantaGenomeCenter', 'RGDS=ExomeSeq', 'RGDT=2016-08-24', 'RGPI=null',
            'RGPG=Test', 'RGPM=Test', 'CREATE_INDEX=true']
    dup_path = '{0}/{1}_MD.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    met_path = '{0}/{1}_MD.metrics'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    mdcmd = ['java', '-Xmx1g', '-jar', picard_path, 'MarkDuplicates', 'I=' + add_path,
             'O=' + dup_path, 'METRICS_FILE=' + met_path, 'REMOVE_DUPLICATES=false',
             'ASSUME_SORTED=true', 'CREATE_INDEX=true']
    fix_path = '{0}/{1}_FM.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    fcmd = ['java', '-Xmx1g', '-jar', picard_path, 'FixMateInformation',
            'I=' + dup_path, 'O=' + fix_path, 'ASSUME_SORTED=true', 'ADD_MATE_CIGAR=true',
            'CREATE_INDEX=true']
    interval_path = '{0}/{1}.intervals'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    trcmd = ['java', '-jar', gatk_path, '-T', 'RealignerTargetCreator', '-o',
             interval_path, '-nt', THREAD_COUNT, '-I', fix_path, '-R', ref_path, '-known',
             dbsnp_path]
    ral_path = '{0}/{1}_IR.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    recmd = ['java', '-jar', gatk_path, '-T', 'IndelRealigner',
             '--targetIntervals', interval_path, '-o', ral_path,
             '-I', fix_path, '-R', ref_path]
    bqs_path = '{0}/{1}.table'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    bqscmd = ['java', '-jar', gatk_path, '-T', 'BaseRecalibrator', '-R', ref_path,
              '-I', ral_path, '-o', bqs_path, '-nct', THREAD_COUNT, '-cov', 'ReadGroupCovariate',
              '-knownSites', dbsnp_path]
    fbam_path = '{0}/{1}_final.bam'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    prcmd = ['java', '-jar', gatk_path, '-T', 'PrintReads', '-R', ref_path, '-I',
             ral_path, '-o', fbam_path, '-BQSR', bqs_path, '-nct', THREAD_COUNT]
    vcf_path = '{0}/variants.vcf'.format(out_path, os.path.splitext(os.path.basename(sam_path))[0])
    hcmd = ['java', '-jar', gatk_path, '-T', 'HaplotypeCaller', '-R', ref_path,
            '-I', fbam_path, '--dbsnp', dbsnp_path, '-o', vcf_path, '-nct', THREAD_COUNT,
            '-gt_mode', 'DISCOVERY']

    #Trim fastq files
    run_process(tcmd, 'Fastq trimming failed')
    #Align the reads using bowtie
    run_process(bcmd, 'Bowtie failed')
    #Add read group information
    run_process(acmd, 'Picard add read groups failed')
    #Mark PCR duplicates
    run_process(mdcmd, 'Picard mark duplicate failed')
    #Fix mate information
    run_process(fcmd, 'Picard fix mate information failed')
    #Run realigner target creator
    run_process(trcmd, 'Realigner Target creator failed')
    #Run indel realigner
    run_process(recmd, 'Indel realigner creator failed')
    #Base quality score recalibration
    run_process(bqscmd, 'Base quality score recalibrator failed')
    #Print Reads
    run_process(prcmd, 'Print reads failed')
    #Haplotype caller
    run_process(hcmd, 'Haplotype caller')

    print('Variant call pipeline completed')
    print('VCF file can be found at {0}'.format(vcf_path))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='VariantCaller')
    parser.add_argument('-t', '--trimmomatic', dest='trim_path', type=str, help='Path to Trimmomatic')
    parser.add_argument('-b', '--bowtie', dest='bowtie_path', type=str, help='Path to Bowtie')
    parser.add_argument('-p', '--picard', dest='picard_path', type=str, help='Path to Picard')
    parser.add_argument('-g', '--gatk', dest='gatk_path', type=str, help='Path to GATK')
    parser.add_argument('-i', '--inputs', dest='input_path', nargs='+', type=str, help='Path to Paired end reads')
    parser.add_argument('-w', '--index', dest='index_path', type=str, help='Path to Reference bowtie index')
    parser.add_argument('-d', '--dbsnp', dest='dbsnp_path', type=str, help='Path to dbSNP vcf file')
    parser.add_argument('-r', '--reference', dest='ref_path', type=str, help='Path to Reference file')
    parser.add_argument('-a', '--adapter', dest='adapter_path', type=str, help='Path to Adapter file')
    parser.add_argument('-o', '--outpath', dest='out_path', type=str, help='Path to Ouput directory')
    args = parser.parse_args()

    main(args.trim_path, args.bowtie_path, args.picard_path, args.gatk_path,
         args.input_path, args.index_path, args.dbsnp_path, args.adapter_path,
         args.ref_path, args.out_path)
