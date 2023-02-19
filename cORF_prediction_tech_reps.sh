#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=500G

SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/orf_prediction_tech_rep
db_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/db/
mkdir -p $output_path

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
    	getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_BPY2.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_BPY2.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
      makeblastdb -in $output_path/$SAMPLE\_BPY2\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_BPY2\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_BPY2\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_BPY2\_query_human_Y_prot_START_STOP_tabular_sorted.txt
     
done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_CDY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_CDY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
      makeblastdb -in $output_path/$SAMPLE\_CDY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_CDY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_CDY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_CDY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_CDY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_CDY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_DAZ.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_DAZ.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
     
      makeblastdb -in $output_path/$SAMPLE\_DAZ\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_DAZ\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_DAZ\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_DAZ\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_HSFY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_HSFY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
      makeblastdb -in $output_path/$SAMPLE\_HSFY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_HSFY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_HSFY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_HSFY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_PRY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_PRY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
      makeblastdb -in $output_path/$SAMPLE\_PRY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_PRY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_PRY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_PRY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_PRY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_PRY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_RBMY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_RBMY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta

      makeblastdb -in $output_path/$SAMPLE\_RBMY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_RBMY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_RBMY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_RBMY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_TSPY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_TSPY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
 
      makeblastdb -in $output_path/$SAMPLE\_TSPY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_TSPY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_TSPY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_TSPY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_VCY_chimp_bonobo.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_VCY_chimp_bonobo.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
 
      makeblastdb -in $output_path/$SAMPLE\_VCY_chimp_bonobo\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_VCY_chimp_bonobo\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_VCY_chimp_bonobo\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_VCY_chimp_bonobo\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_VCY_human.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_VCY_human.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta
      
      makeblastdb -in $output_path/$SAMPLE\_VCY_human\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_VCY_human\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_VCY_human\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_VCY_human\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done

for SAMPLE in 1_1 1_2 2_1 2_2 3_1 3_2 4_1 4_2 5_1 5_2 6_1 6_2 7_1 7_2 8_1 8_2
do
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_XKRY.fa -find 1 -minsize 150 -outseq $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_START_STOP.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_START_STOP.fasta > $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_linearized_START_STOP.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_linearized_START_STOP.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_sorted_START_STOP.fasta
      
      getorf -sequence $SAMPLE_path/$SAMPLE/final_candidates_XKRY.fa -find 3 -minsize 150 -outseq $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_START_STOP_nuc.fasta
      awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_START_STOP_nuc.fasta > $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta
      awk -F '\t' '{printf("%d\t%s\t%s\n",length($2),substr($1,1,15),$0);}' $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_linearized_START_STOP_nuc.fasta | sort -k2,2 -k1,1n | cut -f3,4 | tr "\t" "\n" > $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_sorted_START_STOP_nuc.fasta

      makeblastdb -in $output_path/$SAMPLE\_XKRY\_50aa_and_more_orf_START_STOP.fasta -title pdbaa1 -dbtype prot -out $output_path/$SAMPLE\_XKRY\_pdbaa1 -parse_seqids
      blastp -query $db_path/Ychr_proteins_with_XKRY_mod.fa -db $output_path/$SAMPLE\_XKRY\_pdbaa1 -evalue 1e-30 -outfmt 11 -out $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_results_blastp_START_STOP.txt
      blast_formatter -archive $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_results_blastp_START_STOP.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_tabular.txt
      awk '!x[$2]++' $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_unique.txt
      sort -rn -k4 $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_unique.txt > $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_unique_sorted.txt
      sort -rn -k4 $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_tabular.txt > $output_path/$SAMPLE\_XKRY\_query_human_Y_prot_START_STOP_tabular_sorted.txt

done



