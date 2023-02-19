library(ggplot2)
library(dplyr)
library(magrittr)
library(ggtranscript)
library("RColorBrewer")

#Figure 1A
all_shared_summary_2 %>% filter(!(sample %in% c("Human 2", "Chimpanzee 1"))) -> excluding_2_samples
excluding_2_samples$sample <- factor(excluding_2_samples$sample, levels=c("Human 1", "Chimpanzee 2", "Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))

Fig_1A <- ggplot(excluding_2_samples) +
	aes(x = gene_family, fill = sample) +
	geom_bar(position="dodge", color = “grey”, stat = "count") +
	labs(x = "Gene family", y = "Transcript number") +
	theme_bw() + stat_count(geom = "text", colour = "black", size = 2.5, aes(label = ..count..), position = position_dodge(width = 1), vjust=-1)

ggsave("Fig_1A.png", Fig_1A,height=10,width=10)

#Figure 1B
Fig_1B <- ggplot(excluding_2_samples) +
    aes(x = gene_family, y = transcript_length, fill = sample) +
    geom_boxplot() +
    scale_fill_hue(direction = 1) +
    labs(x = "Gene family", y = "Transcript length") +
    theme_bw() + stat_summary(position = position_dodge(width = 0.75), fun="mean",color="dark grey", size = 0.09)

ggsave("Fig_1B.png", Fig_1B,height=10,width=10)

#Figure 2
All_proteins_50_and_more_sorted_readlength %>% filter(!(sample %in% c("Human 2", "Chimpanzee 1"))) -> Excluding_2_samples_proteins_50_and_more_sorted

Excluding_2_samples_proteins_50_and_more_sorted$sample <- factor(Excluding_2_samples_proteins_50_and_more_sorted$sample, levels=c("Human 1", "Chimpanzee 2", "Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))

Fig_2A <- ggplot(Excluding_2_samples_proteins_50_and_more_sorted) +
    aes(x = orf_length, fill = sample) +
    geom_histogram(bins = 30L) +
    scale_fill_hue(direction = 1) +
    labs(x = "ORF length", y = "ORF count") +
    theme_minimal() +
    facet_wrap(vars(gene_family), scales = "fixed", nrow = 5L)

ggsave("Fig_2A.png", Fig_2A,height=10,width=10)

#Figure 2B
Excluding_2_samples_proteins_50_and_more_sorted$sample <- factor(Excluding_2_samples_proteins_50_and_more_sorted$sample, levels=c("Human 1", "Chimpanzee 2", "Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))

Fig_2B <- ggplot(Excluding_2_samples_proteins_50_and_more_sorted) +
  	aes(x = gene_family, fill = gene_family) +
  	geom_bar(stat = "count") +
  	scale_fill_hue(direction = 1) +
  	labs(x = "Gene family", y = "ORF count") +
  	theme_bw() + stat_count(geom = "text", colour = "black", size = 2.5, aes(label = ..count..), vjust=-0.5) +
  	facet_wrap(vars(sample))

ggsave("Fig_2B.png", Fig_2B,height=10,width=10)

#Figure 3
Human_Y_protein_homologs$sample <- factor(Human_Y_protein_homologs$sample, levels=c("Human 1", "Chimpanzee 2", "Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))

Fig_3 <- ggplot(Human_Y_protein_homologs) +
    aes(x = sample, fill = sample) +
    geom_bar(stat = "count") +
    labs(x = "Great ape species", y = "Number of ORF homologs") +
	theme_bw() + stat_count(geom = "text", colour = "black", size = 2.5, aes(label = ..count..), vjust=-0.5) + facet_wrap(vars(gene_family), nrow = 3L) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylim(0, 50)

ggsave("Fig_3.png", Fig_3,height=10,width=10)

#Figure 4
Fig_4A <- ggplot(Human_Y_protein_homologs) +
  aes(x = alignment_length, fill = sample) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(direction = 1) +
  labs(
    x = "Alignment length (aa) of ORF homologous to human Y proteins",
    y = "Number of homologous to human Y proteins"
  ) +
  theme_bw() +
  facet_wrap(vars(gene_family), scales = "fixed”) 

ggsave("Fig_4A.png", Fig_4A,height=10,width=10)

Fig_4B <- ggplot(Human_Y_protein_homologs) +
  aes(
    x = query_coverage,
    y = sequence_identity,
    colour = sample
  ) +
  geom_point(shape = "circle", size = 1) +
  scale_color_hue(direction = 1) +
  labs(
    x = "Coverage of homologous to human Y proteins",
    y = "Sequence identity of homologous to human Y proteins"
  ) +
  theme_bw() +
  facet_wrap(vars(gene_family), scales = “fixed”)

ggsave("Fig_4B.png", Fig_4B,height=10,width=10)


#Figures 5-14 Alternative splicing patterns of YAG per gene family across great apes.

BPY_All_exons <- BPY2_transcripts_bed%>% dplyr::filter(type == "exon")
BPY_All_exons$species <- factor(BPY_All_exons$species, levels=c("Human reference", "Human 1", "Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
BPY_All_exons$transcript_name <- factor(BPY_All_exons$transcript_name, levels = rev(unique(BPY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("BPY2_CCDS14800.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43")

BPY2_All <- BPY_All_exons %>%
    ggplot(aes(
        xstart = start,
         xend = end,
       y = transcript_name
    )) + scale_y_discrete(labels = rev(my_labels)) +
    geom_range(
        aes(fill = species)
    ) + labs(x = "hg38Y human BPY2 gene coordinates", y = "Transcript name") +
    geom_intron(
        data = to_intron(BPY_All_exons, "transcript_name"),
        aes(strand = strand) 
    )

ggsave("BPY2_All.png",BPY2_All,height=7,width=15)

CDY_All_exons <- CDY_transcripts_bed%>% dplyr::filter(type == "exon")
CDY_All_exons$species <- factor(CDY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
CDY_All_exons$transcript_name <- factor(CDY_All_exons$transcript_name, levels = rev(unique(CDY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("CDY1A_CCDS14801.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56")

CDY_All <- CDY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human CDY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(CDY_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("CDY_All.png",CDY_All,height=7,width=15)

DAZ_All_exons <- DAZ_transcripts_bed%>% dplyr::filter(type == "exon")
DAZ_All_exons$species <- factor(DAZ_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
DAZ_All_exons$transcript_name <- factor(DAZ_All_exons$transcript_name, levels = rev(unique(DAZ_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("DAZ1_CCDS48209.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81", "transcript_82", "transcript_83", "transcript_84", "transcript_85", "transcript_86", "transcript_87", "transcript_88", "transcript_89", "transcript_90", "transcript_91", "transcript_92", "transcript_93", "transcript_94", "transcript_95", "transcript_96", "transcript_97", "transcript_98", "transcript_99", "transcript_100", "transcript_101", "transcript_102", "transcript_103", "transcript_104", "transcript_105", "transcript_106", "transcript_107", "transcript_108", "transcript_109", "transcript_110", "transcript_111", "transcript_112", "transcript_113", "transcript_114", "transcript_115", "transcript_116", "transcript_117", "transcript_118", "transcript_119", "transcript_120", "transcript_121", "transcript_122", "transcript_123","transcript_124", "transcript_125", "transcript_126", "transcript_127", "transcript_128", "transcript_129", "transcript_130","transcript_131","transcript_132","transcript_133","transcript_134","transcript_135","transcript_136","transcript_137","transcript_138","transcript_139","transcript_140","transcript_141","transcript_142","transcript_143", "transcript_144", "transcript_145", "transcript_146", "transcript_147", "transcript_148", "transcript_149", "transcript_150", "transcript_151", "transcript_152", "transcript_153", "transcript_154", "transcript_155", "transcript_156", "transcript_157", "transcript_158", "transcript_159", "transcript_160", "transcript_161", "transcript_162")

DAZ_All <- DAZ_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human DAZ gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(DAZ_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("DAZ_All.png",DAZ_All,height=20,width=20)

HSFY_All_exons <- HSFY_transcripts_bed%>% dplyr::filter(type == "exon")
HSFY_All_exons$species <- factor(HSFY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
HSFY_All_exons$transcript_name <- factor(HSFY_All_exons$transcript_name, levels = rev(unique(HSFY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("HSFY1_CCDS35475.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81", "transcript_82", "transcript_83", "transcript_84", "transcript_85", "transcript_86", "transcript_87", "transcript_88", "transcript_89", "transcript_90", "transcript_91", "transcript_92", "transcript_93", "transcript_94", "transcript_95", "transcript_96", "transcript_97", "transcript_98", "transcript_99", "transcript_100", "transcript_101", "transcript_102", "transcript_103", "transcript_104", "transcript_105", "transcript_106", "transcript_107", "transcript_108", "transcript_109")

HSFY_All <- HSFY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human HSFY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(HSFY_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("HSFY_All.png",HSFY_All,height=15,width=15)

PRY_All_exons <- PRY_transcripts_bed%>% dplyr::filter(type == "exon")
PRY_All_exons$species <- factor(PRY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
PRY_All_exons$transcript_name <- factor(PRY_All_exons$transcript_name, levels = rev(unique(PRY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("PRY_CCDS14799.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7”)

PRY_All <- PRY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human PRY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(PRY_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("PRY_All.png",PRY_All,height=4,width=10)

RBMY_All_exons <- RBMY_transcripts_bed%>% dplyr::filter(type == "exon")
RBMY_All_exons$species <- factor(RBMY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
RBMY_All_exons$transcript_name <- factor(RBMY_All_exons$transcript_name, levels = rev(unique(RBMY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("RBMY1A_CCDS14796.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81", "transcript_82", "transcript_83", "transcript_84", "transcript_85", "transcript_86", "transcript_87", "transcript_88", "transcript_89", "transcript_90", "transcript_91", "transcript_92", "transcript_93", "transcript_94", "transcript_95", "transcript_96", "transcript_97", "transcript_98", "transcript_99", "transcript_100", "transcript_101", "transcript_102", "transcript_103", "transcript_104", "transcript_105", "transcript_106", "transcript_107", "transcript_108", "transcript_109", "transcript_110", "transcript_111", "transcript_112", "transcript_113", "transcript_114", "transcript_115", "transcript_116", "transcript_117", "transcript_118", "transcript_119", "transcript_120", "transcript_121", "transcript_122", "transcript_123","transcript_124", "transcript_125", "transcript_126", "transcript_127", "transcript_128", "transcript_129", "transcript_130","transcript_131","transcript_132","transcript_133","transcript_134","transcript_135","transcript_136","transcript_137","transcript_138","transcript_139","transcript_140","transcript_141","transcript_142","transcript_143", "transcript_144", "transcript_145")

RBMY_All <- RBMY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human RBMY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(RBMY_All_exons, "transcript_name"),
         aes(strand = strand) 
   )

ggsave("RBMY_All.png",RBMY_All,height=20,width=15)

TSPY_All_exons <- TSPY_transcripts_bed%>% dplyr::filter(type == "exon")
TSPY_All_exons$species <- factor(TSPY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
TSPY_All_exons$transcript_name <- factor(TSPY_All_exons$transcript_name, levels = rev(unique(TSPY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("TSPY1_CCDS76071.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81", "transcript_82", "transcript_83", "transcript_84", "transcript_85", "transcript_86", "transcript_87", "transcript_88", "transcript_89", "transcript_90", "transcript_91", "transcript_92", "transcript_93", "transcript_94", "transcript_95", "transcript_96", "transcript_97", "transcript_98", "transcript_99", "transcript_100", "transcript_101", "transcript_102", "transcript_103", "transcript_104", "transcript_105", "transcript_106", "transcript_107", "transcript_108", "transcript_109", "transcript_110", "transcript_111", "transcript_112", "transcript_113", "transcript_114", "transcript_115", "transcript_116", "transcript_117", "transcript_118", "transcript_119", "transcript_120", "transcript_121")

TSPY_All <- TSPY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human TSPY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(TSPY_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("TSPY_All.png",TSPY_All,height=15,width=15)

VCY_h_All_exons <- VCY_h_transcripts_bed%>% dplyr::filter(type == "exon")
VCY_h_All_exons$species <- factor(VCY_h_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
VCY_h_All_exons$transcript_name <- factor(VCY_h_All_exons$transcript_name, levels = rev(unique(VCY_h_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("VCY_CCDS56617.1", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81")

VCY_h_All <- VCY_h_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human VCY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(VCY_h_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("VCY_h_All.png",VCY_h_All,height=15,width=15)

VCY_cb_All_exons <- VCY_cb_transcripts_bed%>% dplyr::filter(type == "exon")
VCY_cb_All_exons$species <- factor(VCY_cb_All_exons$species, levels=c("Chimpanzee reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
VCY_cb_All_exons$transcript_name <- factor(VCY_cb_All_exons$transcript_name, levels = rev(unique(VCY_cb_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("VCY_Ptr_ENSPTRT00000055080", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63")

VCY_cb_All <- VCY_cb_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "Pan_tro_3.0 chimpanzee VCY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(VCY_cb_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("VCY_cb_All.png",VCY_cb_All,height=15,width=15)


XKRY_All_exons <- XKRY_transcripts_bed%>% dplyr::filter(type == "exon")
XKRY_All_exons$species <- factor(XKRY_All_exons$species, levels=c("Human reference", "Human 1", "Chimpanzee 2","Bonobo", "Gorilla", "Bornean orangutan", "Sumatran orangutan"))
XKRY_All_exons$transcript_name <- factor(XKRY_All_exons$transcript_name, levels = rev(unique(XKRY_All_exons$transcript_name)), ordered=TRUE)

my_labels <- c("XKRY_NM_004677.2", "transcript_1", "transcript_2", "transcript_3", "transcript_4", "transcript_5", "transcript_6", "transcript_7", "transcript_8", "transcript_9", "transcript_10", "transcript_11", "transcript_12", "transcript_13", "transcript_14", "transcript_15", "transcript_16", "transcript_17", "transcript_18", "transcript_19", "transcript_20","transcript_21", "transcript_22", "transcript_23", "transcript_24", "transcript_25", "transcript_26", "transcript_27", "transcript_28", "transcript_29", "transcript_30", "transcript_31", "transcript_32", "transcript_33", "transcript_34", "transcript_35", "transcript_36", "transcript_37", "transcript_38", "transcript_39", "transcript_40", "transcript_41", "transcript_42", "transcript_43", "transcript_44", "transcript_45", "transcript_46", "transcript_47", "transcript_48", "transcript_49", "transcript_50", "transcript_51", "transcript_52", "transcript_53", "transcript_54", "transcript_55", "transcript_56", "transcript_57", "transcript_58", "transcript_59", "transcript_60", "transcript_61", "transcript_62", "transcript_63", "transcript_64", "transcript_65", "transcript_66", "transcript_67", "transcript_68", "transcript_69", "transcript_70", "transcript_71", "transcript_72", "transcript_73", "transcript_74", "transcript_75", "transcript_76", "transcript_77", "transcript_78", "transcript_79", "transcript_80", "transcript_81", "transcript_82", "transcript_83", "transcript_84", "transcript_85", "transcript_86", "transcript_87", "transcript_88", "transcript_89", "transcript_90", "transcript_91", "transcript_92", "transcript_93", "transcript_94", "transcript_95", "transcript_96", "transcript_97", "transcript_98", "transcript_99", "transcript_100", "transcript_101", "transcript_102”)

XKRY_All <- XKRY_All_exons %>%
     ggplot(aes(
         xstart = start,
         xend = end,
         y = transcript_name
     )) + scale_y_discrete(labels = rev(my_labels)) +
     geom_range(
         aes(fill = species)
     ) + labs(x = "hg38Y human XKRY gene coordinates", y = "Transcript name") +
     geom_intron(
         data = to_intron(XKRY_All_exons, "transcript_name"),
         aes(strand = strand) 
    )

ggsave("XKRY_All.png",XKRY_All,height=15,width=15)



