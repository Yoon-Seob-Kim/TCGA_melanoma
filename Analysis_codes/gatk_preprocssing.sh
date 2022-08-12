#!/bin/bash
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
for i in SRX4042231
do
$BWA/bwa mem -t 4 $DB/human_g1k_v37.fasta $i\_1.fastq $i\_2.fastq > $i\.sam
rm $i\.fastq $i\_1.fastq $i\_2.fastq
$GATK/gatk --java-options "-Xmx4g" SortSam -I=$i\.sam -O=$i\.bam -SO=coordinate --VALIDATION_STRINGENCY=SILENT --CREATE_INDEX=true
rm $i\.sam
$GATK/gatk --java-options "-Xmx4g" MarkDuplicates -I=$i\.bam -O=$i\_dup.bam -M=$i\_dup.met --VALIDATION_STRINGENCY=SILENT --REMOVE_DUPLICATES=true
rm $i\.bam rm $i\.bai
$GATK/gatk --java-options "-Xmx4g" AddOrReplaceReadGroups -I=$i\_dup.bam -O=$i\_dupRemoved.bam -SO=coordinate --VALIDATION_STRINGENCY=SILENT -ID=$i -LB=$i -PL=Illumina -PU=$i -SM=$i --CREATE_INDEX TRUE
rm $i\_dup.bam $i\_dup.met
$GATK/gatk --java-options "-Xmx4g" BaseRecalibrator -R $DB/human_g1k_v37.fasta -I $i\_dupRemoved.bam --known-sites $DB/dbsnp_138.b37.excluding_sites_after_129.vcf --known-sites $DB/1000G_phase1.indels.b37.vcf -O $i\_recal_data.table
$GATK/gatk --java-options "-Xmx4g" ApplyBQSR -R $DB/human_g1k_v37.fasta -I $i\_dupRemoved.bam -bqsr $i\_recal_data.table -O $BAM_DIR/$i\_b37.bam
rm $i\_dupRemoved.bam $i\_dupRemoved.bai $i\_recal_data.table
$GATK/gatk --java-options "-Xmx4g" GetPileupSummaries -R $DB/human_g1k_v37.fasta -I $BAM_DIR/$i\_b37.bam -O $MT_DIR/$i\_pileups.table -V $DB/af-only-gnomad.raw.sites.b37.vcf.gz --intervals $DB/b37_wgs_calling_regions.v1.list
done

