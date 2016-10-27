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

Some rough notes:  
vcf we generated  
many gold vcf's  
project.NIST.hc.snps.indels.vcf  
get ranges from bed  
then grab snps from vcf that are in region  
or use bedtools intersect  
look for specificity/sensitivity  

ID=AD  
ID=QD

GATK variant quality learning and filtering  
https://software.broadinstitute.org/gatk/guide/article?id=1259

```{sh}
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/project.NIST.hc.snps.indels.vcf
brew install bcftools
bgzip project.NIST.hc.snps.indels.vcf
bgzip variants.vcf
tabix -p vcf project.NIST.hc.snps.indels.vcf
tabix -p vcf variants.vcf
bcftools stats project.NIST.hc.snps.indels.vcf variants.vcf > variants.vchk
apt install python-numpy python-scipy python-matplotlib
apt install texlive-binaries texlive-latex-base texlive-latex-recommended texlive-latex-extra
#pip install numpy
#pip install scipy
#pip install matplotlib
#brew install bcftools
#brew cask install mactex
plot-vcftools -p output_dir --no-PDF variants.vchk
```

# best way to split large genome bam file

```{sh}
bamtools split -in file.bam -reference
```

# snippets

```{sh}
cat ahcg_pipeline/breast_ovarian_cancer_genelist.txt | cut -f2 | tail -n26 | cut -d. -f1 | xargs -n1 -I%%% grep %%% hg19_refGene.txt | cut -f3 | sort -u >bc_oc_chroms.txt
```

# VCF quality trimming using VQSR

http://gatkforums.broadinstitute.org/firecloud/discussion/39/variant-quality-score-recalibration-vqsr  
https://software.broadinstitute.org/gatk/guide/article?id=2805  
https://software.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php  

Edit each script and change any paths needed.  
then, generate the recalibration files with first script (15-30min)  
then, trim vcf using second script (1min)  

```{sh}
bash bcoc_vqsr.sh
bash bcoc_apply_recal.sh
```

# matching variants from vcf with clinical risks

```{sh}
python3 compare_clin_with_vcf.py final_variants_trimmed_and_interected.vcf BRCA1_brca_exchange_variants.csv BRCA2_brca_exchange_variants.csv \
| tee brca_clinical_xref.txt

grep -vi benign brca_clinical_xref.txt > brca_clinical_nonbenign_xref.txt
python3 convert_clintxt_to_vcf.py brca_clinical_nonbenign_xref.txt brca_clinical_nonbenign_xref.vcf
```

# coverage calculator

Min generated brca1.join\_final.bed on Vannberg's server.  
he checked it in here: https://github.com/mjaeyi/ahcg\_pipeline/blob/master/brca1.coverage1.bed  
NOTE: the filename is different.

```{sh}
grep 'NM_007298' bcoc_padded.bed > brca1.bed
samtools view -L brca1.bed data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > new.bam
bedtools genomecov -ibam new.bam -bga > na12878.bga.bed
bedtools intersect -loj -a brca1.bed -b na12878.bga.bed -bed > brca1.join_final.bed
awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$8,$8,$4,$10,$6)}' brca1.join_final.bed \
| sed -E -e 's/^chr//' > brca1.final.bed

bedtools intersect -a brca1.final.bed -b brca_clinical_nonbenign_xref.vcf -wo > brca_clinical_nonbenign_final.bed

awk '{printf("%s\t%s\t%s\t%s\t%s\t%s\n",$1,$2,$3,$4,$10,$6)}' brca1.join_final.bed > brca1.depths.bed
python cov.py brca1.depths.bed brca_depth.txt
Rscript draw_depth.R brca_depth.txt brca_depth.png
```

# report generation

http://www.methods.co.nz/asciidoc/userguide.html  
http://powerman.name/doc/asciidoc  

```{sh}
brew install asciidoc
export XML_CATALOG_FILES="/usr/local/etc/xml/catalog"
pip install --ignore-installed --install-option="--install-scripts=/usr/local/bin" dblatex

awk 'BEGIN {FS="\t"} {
gsub(/PG=/,"",$14)
gsub(/\|/,";",$14)
if ($7 !~ /^chr/) $7 = "chr" $7
printf("%s\t%s\t%s\t%s\t%s\n",$4,$5,$7,$8,$14)
}' brca_clinical_nonbenign_final.bed > brca_final_from_bed.txt

grep -i benign brca_final_from_bed.txt > benign_report_body.txt
test -s benign_report_body.txt || echo No Matches Found > benign_report_body.txt
cat benign_report_body.txt | tr '\t' , > benign_report_body.csv

grep -iv benign brca_final_from_bed.txt > nonbenign_report_body.txt
test -s nonbenign_report_body.txt || echo No Matches Found > nonbenign_report_body.txt
cat nonbenign_report_body.txt | tr '\t' , > nonbenign_report_body.csv

cat article-title.txt \
benign-header.txt \
table-header.txt \
benign_report_body.csv \
table-footer.txt \
nonbenign-header.txt \
table-header.txt \
nonbenign_report_body.csv \
table-footer.txt \
depth-chart.txt \
> NA12878-brca1-final-report.txt

a2x -f pdf -d article \
--dblatex-opts '-P doc.layout="mainmatter" -P doc.publisher.show=0' \
NA12878-brca1-final-report.txt
```
