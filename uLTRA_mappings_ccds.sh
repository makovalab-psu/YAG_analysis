#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50G


SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/ultra_cds_probes
mkdir -p $output_path

/nfs/thumper.galaxyproject.org/home/marta/conda/envs/ultra/bin/uLTRA pipeline /nfs/brubeck.bx.psu.edu/scratch6/marta/Hsa_amp_gene_plus_VCY_Ptr.fa /nfs/brubeck.bx.psu.edu/scratch6/marta/Whole_genome_Hsa/hg38_chrY.vcf /nfs/brubeck.bx.psu.edu/scratch6/marta/Hsa_amp_ccds.fa $output_path/CDS_ultra --t 40
samtools view -Su $output_path/CDS_ultra/reads.sam > $output_path/CDS_ultra\_temp_wgs_hsa.bam
samtools sort $$output_path/CDS_ultra\_temp_wgs_hsa.bam > $output_path/CDS_ultra\_wgs_hsa_sorted.bam
samtools index $output_path/CDS_ultra\_wgs_hsa_sorted.bam
/galaxy/home/marta/conda/envs/bedtools/bedtools2/bin/bamToBed -i $$output_path/CDS_ultra\_wgs_hsa_sorted.bam -split > $output_path/CDS_ultra\_split.bed