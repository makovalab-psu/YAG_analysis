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
    for MODEL in null_model alt_model_1 alt_model_2 alt_model_3 alt_model_4 alt_model_5 alt_model_6 alt_model_7 alt_model_8 alt_model_9 alt_model_10 alt_model_11 alt_model_12 alt_model_13 alt_model_14 alt_model_15 alt_model_16 alt_model_17 alt_model_18 alt_model_19 alt_model_20 alt_model_21 alt_model_22 alt_model_23 alt_model_24 alt_model_25 alt_model_26 alt_model_27 alt_model_28 alt_model_29 alt_model_30 alt_model_31 alt_model_32 alt_model_33 alt_model_34 alt_model_35 alt_model_36 alt_model_37 alt_model_38 alt_model_39 alt_model_40 alt_model_41 alt_model_42 alt_model_43 alt_model_44 alt_model_45 alt_model_46 alt_model_47 alt_model_48 alt_model_49 alt_model_50 alt_model_51 alt_model_52 alt_model_53 alt_model_54 alt_model_55 alt_model_56 alt_model_57 alt_model_58 alt_model_59 alt_model_60 alt_model_61 alt_model_62 alt_model_63 alt_model_64 alt_model_65 alt_model_66 alt_model_67 alt_model_68 alt_model_69 alt_model_70 alt_model_71 alt_model_72 alt_model_73 alt_model_74 alt_model_75 alt_model_76 alt_model_77 alt_model_78 alt_model_79 alt_model_80 alt_model_81 alt_model_82 alt_model_83 alt_model_84 alt_model_85 alt_model_86 alt_model_87 alt_model_88 alt_model_89 alt_model_90 alt_model_91 alt_model_92 alt_model_93 alt_model_94 alt_model_95 alt_model_96 alt_model_97 alt_model_98 alt_model_99 alt_model_100 alt_model_101 alt_model_102 alt_model_103 alt_model_104 alt_model_105 alt_model_106 alt_model_107 alt_model_108 alt_model_109 alt_model_110 alt_model_111 alt_model_112 alt_model_113
    do
    	cd $AMP_GENE_path/$GENE/$MODEL/
		/nfs/brubeck.bx.psu.edu/scratch5/marta/paml4.8/bin/codeml *.ctl

    done
done




