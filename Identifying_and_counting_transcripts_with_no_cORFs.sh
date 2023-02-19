#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=500G


SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/orf_prediction_new


for SAMPLE in 1 2 3 4 5 6 7 8
do
    for GENE in BPY2 CDY DAZ HSFY PRY RBMY TSPY VCY_chimp_bonobo VCY_human XKRY
    do
    	grep "^>" $SAMPLE_path/$SAMPLE/$GENE/shared_final_transcripts.fa > $SAMPLE_path/$SAMPLE/$GENE\_$GENE\_transcript_list.txt
    	sed 's|_[^_]*$||' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP_nuc_shortened.fasta
    	grep "^>" $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP_nuc_shortened.fasta > $output_path/$SAMPLE\_$GENE\_orf_list.txt
		diff <(cat $SAMPLE_path/$SAMPLE/$GENE\_$GENE\_transcript_list.txt | grep ">" | sort)  <(cat $output_path/$SAMPLE\_$GENE\_orf_list.txt | grep ">" | sort) | grep "^<" | awk -F\> '{print $2}' > $output_path/$SAMPLE\_$GENE\_unique_list.txt
		cat $output_path/$SAMPLE\_$GENE\_unique_list.txt | wc -l > $output_path/$SAMPLE\_$GENE\_unique_count.txt
		cat $output_path/1_*unique_count.txt > 1_ALL_unique_counts.txt
		cat $output_path/2_*unique_count.txt > 2_ALL_unique_counts.txt
		cat $output_path/3_*unique_count.txt > 3_ALL_unique_counts.txt
		cat $output_path/4_*unique_count.txt > 4_ALL_unique_counts.txt
		cat $output_path/5_*unique_count.txt > 5_ALL_unique_counts.txt
		cat $output_path/6_*unique_count.txt > 6_ALL_unique_counts.txt
		cat $output_path/7_*unique_count.txt > 7_ALL_unique_counts.txt
		cat $output_path/8_*unique_count.txt > 8_ALL_unique_counts.txt
	done
done
