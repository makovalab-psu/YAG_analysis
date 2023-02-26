#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --mem=500G


AMP_GENE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/Selection_tests/Ampliconic_genes/With_macaque

for GENE in DAZ
do
    for MODEL in alt_model_114 /nfs/brubeck.bx.psu.edu/scratch5/marta/paml4.8/bin/codeml *.ctl -t 64
done



