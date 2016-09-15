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

Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html).  
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

Example shell script to run pipeline:
```{sh}
#!/bin/bash
BASE_DIR="/Users/wmartin45/school/biol8803f"
TRIM=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/trimmomatic-0.36.jar
ADAPTER=$BASE_DIR/ahcg_pipeline/lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa
BOWTIE=$BASE_DIR/ahcg_pipeline/lib/bowtie2-2.2.9/bowtie2
PICARD=$BASE_DIR/ahcg_pipeline/lib/picard.jar
GATK=$BASE_DIR/ahcg_pipeline/lib/GenomeAnalysisTK.jar

python3 $BASE_DIR/ahcg_pipeline/ahcg_pipeline.py \
    -t $TRIM \
    -b $BOWTIE \
    -p $PICARD \
    -g $GATK \
    -i $BASE_DIR/test_r1.fastq.gz $BASE_DIR/test_r2.fastq.gz \
    -w $BASE_DIR/resources/genome/hg19 \
    -d $BASE_DIR/resources/dbsnp/dbsnp_138.hg19.vcf.gz \
    -r $BASE_DIR/resources/genome/hg19.fa \
    -a $ADAPTER \
    -o $BASE_DIR/hw1
```

# Genome in a bottle

link acquired from https://github.com/genome-in-a-bottle/giab_latest_release

1. The latest release for NA12878_HG001 is under:

   ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/

2. Fastq files for testing

   http://vannberg.biology.gatech.edu/data/ahcg2016/fq/NA12878_brca_r1.fastq  
   http://vannberg.biology.gatech.edu/data/ahcg2016/fq/NA12878_brca_r2.fastq

4 sets of bam files  

```{sh}
samtools merge out.bam in1.bam in2.bam in3.bam in4.bam
samtools view HG001.hs37d5.300x.bam -L NM_007294.bed -b -o NM_007294.bam
bedtools bamtofastq -i NM_007294.bam -fq  NM_007294.R1.fq -fq2  NM_007294.R2.fq
```

# vcf gold standard comparision

```{sh}
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf
brew install bcftools
bgzip project.NIST.hc.snps.indels.vcf
bgzip variants.vcf
tabix -p vcf project.NIST.hc.snps.indels.vcf
tabix -p vcf variants.vcf
bcftools stats project.NIST.hc.snps.indels.vcf variants.vcf > variants.vchk
#apt install python-numpy python-scipy python-matplotlib
pip install numpy
pip install scipy
pip install matplotlib
plot-vcftools -p output_dir variants.vchk
```
