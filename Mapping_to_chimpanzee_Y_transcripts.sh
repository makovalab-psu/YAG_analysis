#Mapping merged shared transcripts per sample per gene family to publicly available chimpanzee Y chromosome transcripts

#!/bin/bash
#SBATCH --job-name=marta
#SBATCH --output=marta-%j.out
#SBATCH --error=marta-%j.err
#SBATCH -C new
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --mem=500G

makeblastdb -in Y_Pan_tro_transcripts.fa -title chimp_transcripts -dbtype nucl -out chimp_transcripts -parse_seqids

SAMPLE_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters 
output_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/blastn_query_chimp_Y_transcripts
db_path=/nfs/brubeck.bx.psu.edu/scratch6/marta/Capture_great_apes/All_samples/overhang/new_barcodes/IsoCon_15/annotate_clusters/db_chimp_transcripts
mkdir -p $output_path

for SAMPLE in 1 2 3 4 5 6 7 8
do
    for GENE in BPY2 CDY DAZ HSFY PRY RBMY TSPY VCY_chimp_bonobo VCY_human XKRY
    do
    
		blastn -query $db_path/Y_Pan_tro_transcripts.fa -db $SAMPLE_path/$SAMPLE/$GENE/shared_transcripts -out $output_path/$SAMPLE\_$GENE\_reverse_results_blastn.txt -evalue 1e-30 -outfmt 11
		blast_formatter -archive $output_path/$SAMPLE\_$GENE\_reverse_results_blastn.txt -outfmt "6 qseqid sseqid qcovs length pident evalue bitscore mismatch gaps qstart qend sstart send qseq sseq" -out $output_path/$SAMPLE\_$GENE\_results_blastn_query_chimp_Y_transcripts_tabular.txt	
	
  done
done
