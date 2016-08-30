# ahcg_pipeline
Variant calling pipeline for genomic data analysis

## Requirements

1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)

## OSX/homebrew install instructions (optional)

```{sh}
brew install homebrew/science/bowtie2
brew install homebrew/science/trimmomatic
brew install homebrew/science/picard-tools
```

## Reference genome

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)
Or, you can get a reduced version from [prism](www.prism.gatech.edu/~sravishankar9/resources.tar.gz)

## Test data

Use the following protocol to download and prepare test dataset from NIST sample NA12878

```{sh}
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
gzcat NIST7035_TAAGGCGA_L001_R1_001.fastq.gz \
| head -n 100000 | gzip > test_r1.fastq.gz
gzcat NIST7035_TAAGGCGA_L001_R2_001.fastq.gz \
| head -n 100000 | gzip > test_r2.fastq.gz
```

## Index generation

```{sh}
samtools faidx hg19.fa
picard CreateSequenceDictionary R=hg19.fa O=hg19.dict
bowtie2-build hg19.fa hg19
tabix -p vcf dbsnp_138.hg19.vcf.gz
```

## Help

To access help use the following command:

```{sh}
python3 ahcg_pipeline.py -h
```
