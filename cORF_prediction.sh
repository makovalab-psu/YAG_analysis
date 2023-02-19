#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=50G


source activate /nfs/thumper.galaxyproject.org/home/marta/conda/envs/seqkit
SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/orf_prediction_new
db_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/db/
mkdir -p $output_path

conda activate bioawk 

for SAMPLE in 1 2 3 4 5 6 7 8
do
    for GENE in BPY2 CDY DAZ HSFY PRY RBMY TSPY VCY_chimp_bonobo VCY_human XKRY
    do
    	
      getorf -sequence $SAMPLE_path/$SAMPLE/$GENE/shared_final_transcripts.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/$GENE/shared_final_transcripts.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
  	  conda activate bioawk 
  	  bioawk -cfastx '{print $name, length($seq)}' $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_START_STOP.fasta > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_readlength.txt
      cat $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_readlength.txt | cut -f 2 | sort -n | uniq -c > $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_sorted_readlength.txt

    #ORF human Y homologs

      makeblastdb -in $output_path/$SAMPLE\_$GENE\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_$GENE\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_$GENE\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_$GENE\_query_human_Y_prot_START_STOP_tabular_sorted.txt
    
  	done
done
