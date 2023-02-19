#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=500G



SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/orf_prediction_new
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/ultra_transcripts_with_orfs
mkdir -p $output_path

for SAMPLE in 1 2 3 4 5 6 7 8
do
    for GENE in BPY2 CDY DAZ HSFY PRY RBMY TSPY VCY_chimp_bonobo VCY_human XKRY
    do
		/nfs/thumper.galaxyproject.org/home/marta/conda/envs/ultra/bin/uLTRA pipeline /nfs/brubeck.bx.psu.edu/scratch6/marta/Hsa_amp_gene_plus_VCY_Ptr.fa /nfs/brubeck.bx.psu.edu/scratch6/marta/Whole_genome_Hsa/hg38_chrY.vcf $SAMPLE_path/$SAMPLE\_$GENE\_transcripts_with_orfs.fasta $output_path/$SAMPLE\_$GENE\_ultra --t 40
		samtools view -Su $output_path/$SAMPLE\_$GENE\_ultra/indexed.sam > $output_path/$SAMPLE\_$GENE\_ultra\_temp_wgs_hsa.bam
		samtools sort $output_path/$SAMPLE\_$GENE\_ultra\_temp_wgs_hsa.bam > $output_path/$SAMPLE\_$GENE\_ultra\_wgs_hsa_sorted.bam
		samtools index $output_path/$SAMPLE\_$GENE\_ultra\_wgs_hsa_sorted.bam
		/galaxy/home/marta/conda/envs/bedtools/bedtools2/bin/bamToBed -i $output_path/$SAMPLE\_$GENE\_ultra\_wgs_hsa_sorted.bam -split > $output_path/$SAMPLE\_$GENE\_ultra\_split.bed

	done
done
