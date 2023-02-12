"""
    snakemake --keep-going -j 999999 --cluster "sbatch --exclude={cluster.exclude} --mem {cluster.mem} -c {cluster.cpus-per-task} -N {cluster.Nodes}  -t {cluster.runtime} -J {cluster.jobname} --mail-type={cluster.mail_type} --mail-user={cluster.mail}" --cluster-config cluster.json --configfile experiments.json --latency-wait 100 --verbose -n


    1. Gene level clustering: Run isONclust for gene-family level clustering
    2. Transcript level clustering: Run IsoCon with trimmed reads
    3. Cluster to gene family mapping: Get gene family name from IsoCon clusters by applying edlib alignment mapping primers to reads.
    4. Final transcripts: For clusters from the same gene family (inferred in step 5): Get final transcripts supported by both experimental replicates by running ‘isoform_similarity.py’. This is a script I implemented performing (ii) above.


"""

shell.prefix("set -o pipefail; ")
# configfile: "experiments.json"

wildcard_constraints:
    nr_reads="[\d]+",

####################################################
########## standard python functions ###############
####################################################

import re
import os
import errno
import shutil
import glob

def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

rule all:
   input: 
          # expand(config["ROOT_OUT"] + "/annotate_clusters/{sample}_{replicate}/rule_complete.txt", sample = [1,2,3,4,5,6,7,8], replicate = [1,2] ),
          expand(config["ROOT_OUT"] + "/annotate_clusters/{sample}/{gene_fam}/shared_final_transcripts.fa", sample = [1,2,3,4,5,6,7,8], gene_fam = config["GENE_FAM"]),
          summary_file = config["ROOT_OUT"] + "/all_shared_summary.csv"
          # config["ROOT_OUT"] + "/results/BPY/final_pred.fasta", 
          # config["ROOT_OUT"] + "/results/CDY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/DAZ/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/HSFY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/PRY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/RBMY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/TSPY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/VCY/final_pred.fasta",
          # config["ROOT_OUT"] + "/results/XKRY/final_pred.fasta"



rule isONclust:
    input:  fastq = config["ROOT_OUT"] + "/{sample}_{replicate}_flnc.fastq"
    output: time_and_mem =  config["ROOT_OUT"] + "/time_and_mem/{sample}_{replicate}/isONclust_time_and_mem.txt",
            clusters = config["ROOT_OUT"] + "/isONclust/{sample}_{replicate}/final_clusters.tsv" 
    run:
        time = config["GNUTIME"]
        mkdir_p(config["ROOT_OUT"] + "/time_and_mem/{0}_{1}/".format(wildcards.sample, wildcards.replicate) )
        outfolder = config["ROOT_OUT"] + "/isONclust/{0}_{1}/".format(wildcards.sample, wildcards.replicate)
        mkdir_p(outfolder)
        shell("{time} python /galaxy/home/ksahlin/prefix/source/isONclust/isONclust --t 1 --k 11 --w 15 --fastq {input.fastq}  --outfolder {outfolder}  2>&1 | tee {output.time_and_mem}")
    

rule clusters_to_fastq:
    input: fastq = rules.isONclust.input.fastq,
            clusters = rules.isONclust.output.clusters
    output: flag = config["ROOT_OUT"] + "/isONclust/{sample}_{replicate}/rule_complete.txt"  
    run:
        time = config["GNUTIME"]
        outfolder = config["ROOT_OUT"] + "/isONclust/{0}_{1}/fastq/".format(wildcards.sample, wildcards.replicate)
        shell("{time} python /galaxy/home/ksahlin/prefix/source/isONclust/isONclust write_fastq --clusters {input.clusters} --fastq {input.fastq} --outfolder {outfolder} --N 10")
        shell("touch {output.flag}")


rule IsoCon:
    input:  rules.clusters_to_fastq.output.flag
    output:  flag = config["ROOT_OUT"] + "/IsoCon/{sample}_{replicate}/rule_complete.txt" 
    run: 
        time = config["GNUTIME"]
        infolder =  config["ROOT_OUT"] + "/isONclust/{0}_{1}/fastq/".format(wildcards.sample, wildcards.replicate) 
        outfolder = config["ROOT_OUT"] + "/IsoCon/{0}_{1}/".format(wildcards.sample, wildcards.replicate)   
        time_and_mem = config["ROOT_OUT"] + "/time_and_mem/{0}_{1}/IsoCon_time_and_mem.txt".format(wildcards.sample, wildcards.replicate)
        # Usage: run_isocon <isonclust_cluster_fastq_folder> <outfolder> 
        shell("{time} /galaxy/home/ksahlin/prefix/source/amplicon_analysis/./run_isocon {infolder}  {outfolder} 16 2>&1 | tee {time_and_mem}")
        shell("touch {output.flag}")


rule annotate_IsoCon_clusters:
    input: rules.IsoCon.output.flag
    output: transcripts = expand(config["ROOT_OUT"] + "/annotate_clusters/{{sample}}_{{replicate}}/final_candidates_{gene_fam}.fa", gene_fam = config["GENE_FAM"]) 
    run:
        infiles = config["ROOT_OUT"] + "/IsoCon/{0}_{1}/*/final_candidates.fa".format(wildcards.sample, wildcards.replicate)
        probes = "/galaxy/home/ksahlin/prefix/source/amplicon_analysis/probe_sequences.fa"
        outfolder = config["ROOT_OUT"] + "/annotate_clusters/{0}_{1}/".format(wildcards.sample, wildcards.replicate)
        shell("rm -r {outfolder}")
        mkdir_p(outfolder)

        # create empty file to assure all output files are created
        for t in output.transcripts:
            shell("touch {t}")

        # add actual predictions to files
        for f in glob.glob(infiles):
            print(f)
            base=os.path.basename(f)
            prefix = os.path.splitext(base)[0]
            shell("python /galaxy/home/ksahlin/prefix/source/amplicon_analysis/annotate_clusters.py  {f} {probes} {outfolder} {prefix}")
        # shell("touch {output.flag}")

rule merge_clusters:
    input: replicate1 = config["ROOT_OUT"] + "/annotate_clusters/{sample}_1/final_candidates_{gene_fam}.fa",
            replicate2 = config["ROOT_OUT"] + "/annotate_clusters/{sample}_2/final_candidates_{gene_fam}.fa"
    output: shared_transcripts = config["ROOT_OUT"] + "/annotate_clusters/{sample}/{gene_fam}/shared_final_transcripts.fa",
            info_file = config["ROOT_OUT"] + "/annotate_clusters/{sample}/{gene_fam}/merged_info.txt"
    run:
        outfolder = config["ROOT_OUT"] + "/annotate_clusters/{0}/{1}/".format(wildcards.sample, wildcards.gene_fam)
        shell("python /galaxy/home/ksahlin/prefix/source/amplicon_analysis/isoform_similarity.py {input.replicate1} {input.replicate2} {outfolder} shared_final_transcripts > {output.info_file}")


rule output_support_info:
    input: rules.merge_clusters.output.shared_transcripts
    output: summary = config["ROOT_OUT"] + "/annotate_clusters/{sample}/{gene_fam}/shared_final_transcripts_supporting/summary.csv",
    run:
        workfolder = config["ROOT_OUT"] + "/annotate_clusters/{0}/{1}/shared_final_transcripts_supporting/".format(wildcards.sample, wildcards.gene_fam)

        reads_repl1 = config["ROOT_OUT"] + "/{0}_1_flnc.fastq".format(wildcards.sample)
        reads_repl2 = config["ROOT_OUT"] + "/{0}_2_flnc.fastq".format(wildcards.sample)

        # merge all the cluster info files within a sample and replicate to have all consensus available 
        # /nfs/brubeck.bx.psu.edu/scratch4/ksahlin/amplicon_2021/IsoCon/1_1/0/cluster_info.tsv
        tsv_files_repl1 = config["ROOT_OUT"] + "/IsoCon/{0}_1/*/cluster_info.tsv".format(wildcards.sample)
        merged_tsv_repl1 = config["ROOT_OUT"] + "/annotate_clusters/{0}/{1}/merged_tsv_repl1_tmp.tsv".format(wildcards.sample, wildcards.gene_fam)
        shell("> {merged_tsv_repl1}")
        for f in glob.glob(tsv_files_repl1):
            shell(" cat {f} >> {merged_tsv_repl1}")

        tsv_files_repl2 = config["ROOT_OUT"] + "/IsoCon/{0}_2/*/cluster_info.tsv".format(wildcards.sample)
        merged_tsv_repl2 = config["ROOT_OUT"] + "/annotate_clusters/{0}/{1}/merged_tsv_repl2_tmp.tsv".format(wildcards.sample, wildcards.gene_fam)
        shell("> {merged_tsv_repl2}")
        for f in glob.glob(tsv_files_repl2):
            shell(" cat {f} >> {merged_tsv_repl2}")

        infiles =  config["ROOT_OUT"] + "/annotate_clusters/{0}/{1}/shared_final_transcripts_supporting/*.fa".format(wildcards.sample, wildcards.gene_fam)

        shell("> {output.summary}")
        for f in glob.glob(infiles):
            shell("python /galaxy/home/ksahlin/prefix/source/amplicon_analysis/print_stats.py {f} {reads_repl1} {reads_repl2} {merged_tsv_repl1} {merged_tsv_repl2} {workfolder} >> {output.summary}")
        # python print_stats.py final_fa_file, reads_repl1, reads_repl2, tsv_repl1, tsv_repl2, outfolder


rule summary_file:
    input: summary_files = expand(config["ROOT_OUT"] + "/annotate_clusters/{sample}/{gene_fam}/shared_final_transcripts_supporting/summary.csv", sample = [1,2,3,4,5,6,7,8], gene_fam = config["GENE_FAM"]) 
    output: summary_file = config["ROOT_OUT"] + "/all_shared_summary.csv",
    run:
        shell("> {output.summary_file}")
        s = open(config["ROOT_OUT"] + "/all_shared_summary.csv", "w")
        for f in input.summary_files:
            gene_fam = f.split("/")[-3]
            sample = f.split("/")[-4]
            for line in open(f, "r"):
                l = line.strip()
                s.write("{0},{1},{2}\n".format(sample, gene_fam,l))

        s.close()
