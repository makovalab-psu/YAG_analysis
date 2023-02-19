#Mapping merged shared transcripts per sample per gene family to orangutan Y chromosome protein-coding sequences

#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=500G

makeblastdb -in Pab_Y_genes_cds.fa -title orang_CDS -dbtype nucl -out orang_CDS -parse_seqids

SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/blastn_orang_CDS
db_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes//IsoCon_15/annotate_clusters/db_orang_CDS/orang_CDS
mkdir -p $output_path

for SAMPLE in 1 2 3 4 5 6 7 8
do
    for GENE in BPY2 CDY DAZ HSFY PRY RBMY TSPY VCY_chimp_bonobo VCY_human XKRY
    do
    
		  blastn -query $SAMPLE_path/$SAMPLE/$GENE/shared_final_transcripts.fa -db $db_path -out $output_path/$SAMPLE\_$GENE\_results_blastn.txt -num_alignments 1 -evalue 1e-30 -outfmt 11
		  blast_formatter -archive $output_path/$SAMPLE\_$GENE\_results_blastn.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_$GENE\_results_blastn_orang_CDS_tabular.txt	
	  
    done
done
