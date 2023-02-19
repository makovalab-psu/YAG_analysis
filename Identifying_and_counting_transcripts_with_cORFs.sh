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
    	sed -e "s/^\(.*\)/>\1/" $output_path/$SAMPLE\_$GENE\_unique_list.txt > $output_path/$SAMPLE\_$GENE\_unique_list_appended.txt
		diff <(cat $SAMPLE_path/$SAMPLE/$GENE\_$GENE\_transcript_list.txt | grep ">" | sort)  <(cat $output_path/$SAMPLE\_$GENE\_unique_list_appended.txt | grep ">" | sort) | grep "^<" | awk -F\> '{print $2}' > $output_path/$SAMPLE\_$GENE\_transcripts_with_ORFs.txt
		grep -A 1 -wFf $output_path/$SAMPLE\_$GENE\_transcripts_with_ORFs.txt $SAMPLE_path/$SAMPLE/$GENE/shared_final_transcripts.fa > $output_path/$SAMPLE\_$GENE\_transcripts_with_orfs.fasta
		cat $output_path/$SAMPLE\_$GENE\_transcripts_with_ORFs.txt | wc -l > $output_path/$SAMPLE\_$GENE\_transcripts_with_ORFs_count.txt
		cat $output_path/1_*transcripts_with_ORFs_count.txt > 1_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/2_*transcripts_with_ORFs_count.txt > 2_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/3_*transcripts_with_ORFs_count.txt > 3_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/4_*transcripts_with_ORFs_count.txt > 4_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/5_*transcripts_with_ORFs_count.txt > 5_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/6_*transcripts_with_ORFs_count.txt > 6_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/7_*transcripts_with_ORFs_count.txt > 7_ALL_transcripts_with_ORFs_count.txt
		cat $output_path/8_*transcripts_with_ORFs_count.txt > 8_ALL_transcripts_with_ORFs_count.txt
	done
done
