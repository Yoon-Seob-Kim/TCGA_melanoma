#!/bin/bash
#SBATCH -J vcall
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -o %x.o%j
#SBATCH -e %x.e%j
#SBATCH -p mrcs
#SBATCH -w MRC4
#tool_path
BWA=/home/kysbbubbu/tools/bwa-0.7.17
GATK=/home/kysbbubbu/tools/gatk-4.1.3.0
QUALIMAP=/home/kysbbubbu/tools/qualimap_v2.2.1/qualimap_v2.2.1
DB=/home/kysbbubbu/genomicdb

##OUTPUT
FASTQC_DIR=./01_FastQC
BAM_DIR=./02_BAM
BAMQC_DIR=./03_BAMQC
MT_DIR=./04_MT
ANNO_DIR=./05_Annovar
SEQZ_DIR=./06_SEQZ

list1="SRX1065456 SRX1065476 SRX1065485 SRX1065596 SRX1065634 SRX1065653 SRX1065732 SRX1065735 SRX1065737 SRX1065739"
list2="SRX1065457 SRX1065477 SRX1065486 SRX1065597 SRX1065635 SRX1065654 SRX1065733 SRX1065736 SRX1065738 SRX1065740"
echo $list1 | sed 's/ /\n/g' > /tmp/a.$$
echo $list2 | sed 's/ /\n/g' > /tmp/b.$$
paste /tmp/a.$$ /tmp/b.$$ | while read item1 item2; do
$GATK/gatk --java-options "-Xmx4g" Mutect2 -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$item1\_b37.bam -I $BAM_DIR/$item2\_b37.bam -normal $item1 -tumor $item2 --intervals $DB/b37_wgs_calling_regions.v1.list --germline-resource $DB/af-only-gnomad.raw.sites.b37.vcf.gz -O $MT_DIR/$item2\.vcf.gz
$GATK/gatk --java-options "-Xmx4g" CalculateContamination -I $MT_DIR/$item2\_pileups.table -matched $MT_DIR/$item1\_pileups.table -O $MT_DIR/$item2\_contamination.table
$GATK/gatk --java-options "-Xmx4g" FilterMutectCalls -V $MT_DIR/$item2\.vcf.gz --contamination-table $MT_DIR/$item2\_contamination.table -R $DB/human_g1k_v37.fasta -O $MT_DIR/$item2\_filt.vcf
awk -F '\t' '{if($1== "#CHROM") print; else if($7 == "PASS") print}' $MT_DIR/$item2\_filt.vcf > $ANNO_DIR/$item2\.vcf
done
rm /tmp/a.$$
rm /tmp/b.$$
