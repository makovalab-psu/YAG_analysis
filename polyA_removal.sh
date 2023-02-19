#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=500G

conda activate /galaxy/home/marta/conda/envs/isoseq3

isoseq3 refine --require-polya merged.fl.5p--sample_dT1_3p.bam barcodes_1.fasta 1_1_flnc.bam
bamtools convert -format fastq -in 1_1_flnc.bam > 1_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT2_3p.bam barcodes_2.fasta 1_2_flnc.bam
bamtools convert -format fastq -in 1_2_flnc.bam > 1_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT3_3p.bam barcodes_3.fasta 2_1_flnc.bam
bamtools convert -format fastq -in 2_1_flnc.bam > 2_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT3_3p.bam barcodes_3.fasta 2_1_flnc.bam
bamtools convert -format fastq -in 2_1_flnc.bam > 2_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT4_3p.bam barcodes_4.fasta 2_2_flnc.bam
bamtools convert -format fastq -in 2_2_flnc.bam > 2_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT5_3p.bam barcodes_5.fasta 3_1_flnc.bam
bamtools convert -format fastq -in 3_1_flnc.bam > 3_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT6_3p.bam barcodes_6.fasta 3_2_flnc.bam
bamtools convert -format fastq -in 3_2_flnc.bam > 3_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT7_3p.bam barcodes_7.fasta 4_1_flnc.bam
bamtools convert -format fastq -in 4_1_flnc.bam > 4_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT8_3p.bam barcodes_8.fasta 4_2_flnc.bam
bamtools convert -format fastq -in 4_2_flnc.bam > 4_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT9_3p.bam barcodes_9.fasta 5_1_flnc.bam
bamtools convert -format fastq -in 5_1_flnc.bam > 5_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT10_3p.bam barcodes_10.fasta 5_2_flnc.bam
bamtools convert -format fastq -in 5_2_flnc.bam > 5_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT11_3p.bam barcodes_11.fasta 6_1_flnc.bam
bamtools convert -format fastq -in 6_1_flnc.bam > 6_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT12_3p.bam barcodes_12.fasta 6_2_flnc.bam
bamtools convert -format fastq -in 6_2_flnc.bam > 6_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT13_3p.bam barcodes_13.fasta 7_1_flnc.bam
bamtools convert -format fastq -in 7_1_flnc.bam > 7_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT14_3p.bam barcodes_14.fasta 7_2_flnc.bam
bamtools convert -format fastq -in 7_2_flnc.bam > 7_2_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT15_3p.bam barcodes_15.fasta 8_1_flnc.bam
bamtools convert -format fastq -in 8_1_flnc.bam > 8_1_flnc.fastq

isoseq3 refine --require-polya merged.fl.5p--sample_dT16_3p.bam barcodes_16.fasta 8_2_flnc.bam
bamtools convert -format fastq -in 8_2_flnc.bam > 8_2_flnc.fastq
