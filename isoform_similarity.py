
from __future__ import print_function
import os,sys
import argparse
import re
import math
import numpy as np
from collections import defaultdict
import string
import fractions
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
# from matplotlib_venn import venn2, venn2_circles
import edlib
import errno
import shutil


def edlib_traceback(x, y, mode="NW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k)
    ed = result["editDistance"]
    if ed == -1:
        ed= len(x)

    if task == "path":
        locations =  result["locations"]
        cigar =  result["cigar"]
        return (ed, locations, cigar)
    else:
        return (ed, None, None)

'''
    Below function taken from https://github.com/lh3/readfq/blob/master/readfq.py
'''

def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].split()[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, (''.join(seqs), None) # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, (seq, ''.join(seqs)); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, (seq, None) # yield a fasta record instead
                break

# def read_fasta(fasta_file):
#     fasta_seqs = {}
#     k = 0
#     temp = ''
#     accession = ''
#     for line in fasta_file:
#         if line[0] == '>' and k == 0:
#             accession = line[1:].strip()
#             fasta_seqs[accession] = ''
#             k += 1
#         elif line[0] == '>':
#             yield accession, temp
#             temp = ''
#             accession = line[1:].strip()
#         else:
#             temp += line.strip().upper()
#     if accession:
#         yield accession, temp

def check_inexact(set1, set2):
    for seq1 in set1:
        best_ed = len(seq1)
        best_cig = ""
        best_loc = ""
        for seq2 in set2:
            if math.fabs(len(seq1) - len(seq2)) > best_ed:
                continue

            edit_distance, locations, cigar = edlib_traceback(seq1, seq2, mode="HW", task="path", k=min(100, best_ed))
            if edit_distance < best_ed:
                best_ed = edit_distance
                best_cig = cigar
                best_loc = locations

        if best_ed > 0:
            print(best_ed,best_cig, best_loc)
            if best_ed > 500:
                print(seq1)

def check_exact_substrings_between_sets(dict1, dict2, max_five_prime_offset, max_three_prime_offset):

    l1 = sorted([(seq, acc) for seq, acc in dict1.items()], key=lambda x: len(x[0]))
    l2 = sorted([(seq, acc) for seq, acc in dict2.items()], key=lambda x: len(x[0]), reverse=True)
    total_exact_substrings = 0
    is_validated = {}
    all_supporting_transcripts = defaultdict(list)

    for i, (seq1, acc1) in enumerate(l1):
        # if i % 200 == 0:
        #     print("inferring seq", i)
        is_substring = False 
        is_exact = False 
        for j, (seq2, acc2) in enumerate(l2):

            if seq1 == seq2:
                is_exact = True
                acc_superstring = acc2
                superstring = seq2
                all_supporting_transcripts[superstring].append((acc1, seq1))
                break
            elif seq1 in seq2:
                start_p = seq2.index(seq1)
                len_superstring = len(seq2)
                acc_superstring = acc2
                superstring = seq2
                end_offset = len(seq2) - start_p - len(seq1)
                beg_offset = start_p
                is_substring = True 
                if beg_offset <= max_five_prime_offset and end_offset <= max_three_prime_offset:
                    all_supporting_transcripts[superstring].append((acc1, seq1))
                # print("tmp....substring:", "beginning offset:", start_p, "end offset:", end_offset, "seq1_len:", len(seq1), "seq2_len:", len_superstring, acc1, acc_superstring)

            if len(seq1) > len(seq2):
                break

        if is_exact:
            # print("EXACT MATCH:", acc1, acc2)
            total_exact_substrings += 1
            is_validated[superstring] = acc_superstring

        elif is_substring:
            # print("SUBSTRING:", start_p, len(seq1), len_superstring, acc1, acc_superstring)
            if beg_offset <= max_five_prime_offset and end_offset <= max_three_prime_offset:
                # total_exact_substrings += 1
                is_validated[superstring] = acc_superstring
                # print("Within threshold!")
    # print(len(is_validated), len(all_supporting_transcripts))
    # assert len(is_validated) == len(all_supporting_transcripts)
    return total_exact_substrings, is_validated, all_supporting_transcripts


def check_exact_substrings_within_sets(dict1, max_five_prime_offset, max_three_prime_offset):
    l1 = sorted([(seq, acc) for seq, acc in dict1.items()], key=lambda x: len(x[0]))
    l1_rev = sorted([(seq, acc) for seq, acc in dict1.items()], key=lambda x: len(x[0]), reverse=True)
    total_exact_substrings = 0
    # is_validated = {}
    all_merged_transcripts = defaultdict(list)
    for i, (seq1, acc1) in enumerate(l1):
        # if i % 200 == 0:
        #     print("inferring seq", i)
        merge = False 
        for j, (seq2, acc2) in enumerate(l1_rev):

            if seq1 == seq2:
                break
            elif (seq1 in seq2) and (len(seq1) < len(seq2)):
                start_p = seq2.index(seq1)
                end_offset = len(seq2) - start_p - len(seq1)
                beg_offset = start_p
                if beg_offset <= max_five_prime_offset and end_offset <= max_three_prime_offset:
                    len_superstring = len(seq2)
                    acc_superstring = acc2
                    superstring = seq2
                    merge = True 
                    # print("tmp....substring:", "beginning offset:", start_p, "end offset:", end_offset, "seq1_len:", len(seq1), "seq2_len:", len_superstring, acc1, acc_superstring)

            if len(seq1) > len(seq2):
                break


        if merge:
            total_exact_substrings += 1
            all_merged_transcripts[superstring].append((acc1, acc_superstring, seq1))
            # is_validated[superstring] = acc_superstring
            # print("SUBSTRING", beg_offset, end_offset, acc1, len(seq1),"bp", " merged into :", acc_superstring, len(superstring), "bp." )
        else: # not merged, original sequence, seq1, and points to itself
            all_merged_transcripts[seq1].append((acc1, acc1, seq1))

    return total_exact_substrings, all_merged_transcripts


def merge_transcripts(transcripts):
    non_redundant = defaultdict(list)
    is_merged = {}
    tr = [(seq, l) for seq, l in transcripts.items()]
    for (seq, l) in tr[::-1]:
        for (sub_acc, super_acc, _) in l:
            if super_acc not in is_merged:
                non_redundant[super_acc].append(sub_acc)
                is_merged[sub_acc] = super_acc
            else:
                longest_tr = is_merged[super_acc]
                non_redundant[longest_tr].append(sub_acc)
                is_merged[sub_acc] = longest_tr

    return non_redundant


def mkdir_p(path):
    print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def main(args):
    sample1_seqs = {seq.upper() : acc for acc, (seq, qual) in  readfq(open(args.sample1, 'r'))}
    sample1_acc_to_seqs = {acc: seq.upper() for acc, (seq, qual) in  readfq(open(args.sample1, 'r'))}

    final_supps1 = []
    for a in sample1_seqs.values():
        final_supp = a.split("_")[4] 
        final_supps1.append(int(final_supp))


    sample2_seqs = {seq.upper():acc for acc, (seq, qual) in  readfq(open(args.sample2, 'r'))}
    sample2_acc_to_seqs = {acc: seq.upper() for acc, (seq, qual) in  readfq(open(args.sample2, 'r'))}
    final_supps2 = []
    for a in sample2_seqs.values():
        final_supp = a.split("_")[4] 
        final_supps2.append(int(final_supp))


    print(len(sample1_seqs), len(sample2_seqs))
    # r = venn2([sample1_seqs.keys(), sample2_seqs.keys()], (args.names[0], args.names[1]))
    # plt.savefig(os.path.join(args.outfolder, args.plot_prefix +'_perfect_ss.pdf' ) )

    exact_matches = set(sample1_seqs.keys()) & set(sample2_seqs.keys())
    print("Exact strings between sets:", len(exact_matches) )
    shared_supps1 = []
    shared_supps2 = []
    for seq in exact_matches:
        # print("identical")
        # print(sample1_seqs[seq])
        # print(sample2_seqs[seq])
        # print()
        shared_supp = int(sample1_seqs[seq].split("_")[4])
        shared_supps1.append(shared_supp)
        shared_supp = int(sample2_seqs[seq].split("_")[4])
        shared_supps2.append(shared_supp)

    
    print("All    supports sample 1:", sorted(final_supps1, reverse=True))
    print("Shared supports sample 1:", sorted(shared_supps1, reverse=True))

    print("All    supports sample 2:", sorted(final_supps2, reverse=True))
    print("Shared supports sample 2:", sorted(shared_supps2, reverse=True))

    # fasta_final = open(os.path.join(args.outfolder, "final.fasta"), "w")
    # for seq in sorted(set(sample1_seqs.keys()).union(set(sample2_seqs.keys())), key=len, reverse=True):
    #     if seq in sample1_seqs and seq in sample2_seqs:
    #         acc = "SAMPLE1_" + sample1_seqs[seq] + "SAMPLE2_" + sample2_seqs[seq] 
    #         fasta_final.write(">{0}\n{1}\n".format(acc, seq))
    #     elif seq in sample1_seqs:
    #         supp = int(sample1_seqs[seq].split("_")[4])
    #         if supp > 10:
    #             fasta_final.write(">{0}\n{1}\n".format(sample1_seqs[seq], seq))
    #     elif seq in sample2_seqs:
    #         supp = int(sample2_seqs[seq].split("_")[4])
    #         if supp > 10:
    #             fasta_final.write(">{0}\n{1}\n".format(sample2_seqs[seq], seq))

    # fasta_final.close()


    total_exact_substrings, is_validated_sample2, supporting_seqs_to_sample2 = check_exact_substrings_between_sets(sample1_seqs, sample2_seqs, args.merge_5prime, args.merge_3prime)
    print("Exact substrings of set1 to set2:", total_exact_substrings)
    print("TOTAL SEQS VALIDATED OTHER EXPERIMENTAL SAMPLE:", len(is_validated_sample2))
    print()

    total_exact_substrings, is_validated_sample1, supporting_seqs_to_sample1 = check_exact_substrings_between_sets(sample2_seqs, sample1_seqs, args.merge_5prime, args.merge_3prime)
    print("Exact substrings of set2 to set1:", total_exact_substrings)
    print("TOTAL SEQS VALIDATED OTHER EXPERIMENTAL SAMPLE:", len(is_validated_sample1))
    print()
    exact_matches = set(is_validated_sample2.keys()) & set(is_validated_sample1.keys())
    print("Exact unique strings between validates seqs:", len(exact_matches) )


    # fasta_final = open(os.path.join(args.outfolder, args.prefix + ".fa"), "w")
    # indiv_fa_folder = args.outfolder + "/" + args.prefix + "_supporting"
    # shutil.rmtree(indiv_fa_folder, ignore_errors=True)
    # mkdir_p(indiv_fa_folder)
    # is_validated = 0
    # for seq in sorted(set(is_validated_sample2.keys()).union(set(is_validated_sample1.keys())), key=len, reverse=True):
    #     if seq in is_validated_sample2 and seq in is_validated_sample1:
    #         acc1 = "_".join(is_validated_sample1[seq].split("_")[:4])
    #         acc2 = "_".join(is_validated_sample2[seq].split("_")[:4])
    #         acc = "REPLICATE1_" + acc1 + "REPLICATE2_" + acc2 
    #         fasta_final.write(">{0}\n{1}\n".format(acc, seq))
    #         is_validated += 1
    #         indv_fasta = open(os.path.join(indiv_fa_folder, acc + ".fa"), "w")
    #         indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
    #         for (supp_acc, supp_seq) in supporting_seqs_to_sample2[seq]:
    #             indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + supp_acc, supp_seq))
    #         for (supp_acc, supp_seq) in supporting_seqs_to_sample1[seq]:
    #             indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + supp_acc, supp_seq))
    #         indv_fasta.close()
    #     elif seq in is_validated_sample2:
    #         acc2 = "_".join(is_validated_sample2[seq].split("_")[:4])
    #         acc =  "longer_from_REPLICATE2_" + acc2
    #         fasta_final.write(">{0}\n{1}\n".format(acc, seq))
    #         is_validated += 1
    #         indv_fasta = open(os.path.join(indiv_fa_folder, acc + ".fa"), "w")
    #         indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
    #         for (supp_acc, supp_seq) in supporting_seqs_to_sample2[seq]:
    #             indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + supp_acc, supp_seq))

    #     elif seq in is_validated_sample1:
    #         acc1 = "_".join(is_validated_sample1[seq].split("_")[:4])
    #         acc =  "longer_from_REPLICATE1_" + acc1
    #         fasta_final.write(">{0}\n{1}\n".format(acc, seq))
    #         is_validated += 1
    #         indv_fasta = open(os.path.join(indiv_fa_folder, acc + ".fa"), "w")
    #         indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
    #         for (supp_acc, supp_seq) in supporting_seqs_to_sample1[seq]:
    #             indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + supp_acc, supp_seq))
    # fasta_final.close()

    is_validated_both = {**is_validated_sample1, **is_validated_sample2}
    # print()
    # print()
    _, supported_merged_transcripts = check_exact_substrings_within_sets(is_validated_both, args.merge_5prime, args.merge_3prime)
    final_merged = merge_transcripts(supported_merged_transcripts)
    # for x,y in final_merged.items():
    #     print(x,y)
    is_validated_both_acc_to_seq =  {acc: seq for seq, acc in is_validated_both.items()}
    # create dict of final merged transcripts: seq -> acc
    seq_to_merged = {is_validated_both_acc_to_seq[acc]: acc for acc, l in final_merged.items()}

    fasta_final_m = open(os.path.join(args.outfolder, args.prefix + ".fa"), "w")
    indiv_fa_folder_m = args.outfolder + "/" + args.prefix + "_supporting"
    shutil.rmtree(indiv_fa_folder_m, ignore_errors=True)
    mkdir_p(indiv_fa_folder_m)
    for seq in sorted(seq_to_merged.keys(), key=len, reverse=True):
        has_been_printed = set()
        final_tr_acc = seq_to_merged[seq]
        list_of_merged_tr = final_merged[final_tr_acc]
        if seq in is_validated_sample2 and seq in is_validated_sample1:
            acc1 = "_".join(is_validated_sample1[seq].split("_")[:4])
            acc2 = "_".join(is_validated_sample2[seq].split("_")[:4])
            acc = "REPLICATE1_" + acc1 + "REPLICATE2_" + acc2 
            fasta_final_m.write(">{0}\n{1}\n".format(acc, seq))
            indv_fasta = open(os.path.join(indiv_fa_folder_m, acc + ".fa"), "w")
            indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
            has_been_printed.add(is_validated_sample1[seq])
            has_been_printed.add(is_validated_sample2[seq])

            for (supp_acc, supp_seq) in supporting_seqs_to_sample2[seq]:
                indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + supp_acc, supp_seq))
                has_been_printed.add(supp_acc)
            for (supp_acc, supp_seq) in supporting_seqs_to_sample1[seq]:
                indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + supp_acc, supp_seq))
                has_been_printed.add(supp_acc)
            
            for tr_acc in list_of_merged_tr:
                if tr_acc not in has_been_printed:
                    # print("BOTH HERE!", tr_acc)
                    if tr_acc in sample1_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + tr_acc, sample1_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    elif tr_acc in sample2_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + tr_acc, sample2_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    else:
                        print("BUG!!!", tr_acc)
                else:
                    pass
                    # print("Has been printed")

            indv_fasta.close()

        elif seq in is_validated_sample2:
            acc2 = "_".join(is_validated_sample2[seq].split("_")[:4])
            acc =  "longer_from_REPLICATE2_" + acc2
            fasta_final_m.write(">{0}\n{1}\n".format(acc, seq))
            indv_fasta = open(os.path.join(indiv_fa_folder_m, acc + ".fa"), "w")
            indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
            has_been_printed.add(is_validated_sample2[seq])

            for (supp_acc, supp_seq) in supporting_seqs_to_sample2[seq]:
                indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + supp_acc, supp_seq))
                has_been_printed.add(supp_acc)

            for tr_acc in list_of_merged_tr:
                if tr_acc not in has_been_printed:
                    # print("longer_from_REPLICATE2 HERE!", tr_acc)
                    if tr_acc in sample1_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + tr_acc, sample1_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    elif tr_acc in sample2_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + tr_acc, sample2_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    else:
                        print("BUG!!!", tr_acc)
                else:
                    pass
                    # print("Has been printed")

        elif seq in is_validated_sample1:
            acc1 = "_".join(is_validated_sample1[seq].split("_")[:4])
            acc =  "longer_from_REPLICATE1_" + acc1
            fasta_final_m.write(">{0}\n{1}\n".format(acc, seq))
            indv_fasta = open(os.path.join(indiv_fa_folder_m, acc + ".fa"), "w")
            indv_fasta.write(">{0}\n{1}\n".format(acc, seq))
            has_been_printed.add(is_validated_sample1[seq])

            for (supp_acc, supp_seq) in supporting_seqs_to_sample1[seq]:
                indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + supp_acc, supp_seq))
                has_been_printed.add(supp_acc)

            for tr_acc in list_of_merged_tr:
                if tr_acc not in has_been_printed:
                    # print("longer_from_REPLICATE1 HERE!", tr_acc)
                    if tr_acc in sample1_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE1_" + tr_acc, sample1_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    elif tr_acc in sample2_acc_to_seqs:
                        indv_fasta.write(">{0}\n{1}\n".format("REPLICATE2_" + tr_acc, sample2_acc_to_seqs[tr_acc]))
                        has_been_printed.add(supp_acc)
                    else:
                        print("BUG!!!", tr_acc)
                else:
                    pass
                    # print("Has been printed")
    fasta_final_m.close()

    # what remains is to print the transcripts in final_merged to the separate outfiles instead of the non-merged above!
    # However, each non-merged transcript above in turn, contains supporting merged transcripts... Need to be carefull here so we get the complete
    # set of transcripts

    # print()
    # print()
    _, sample1_merged_transcripts = check_exact_substrings_within_sets(sample1_seqs, args.merge_5prime, args.merge_3prime)
    sample1_merged = merge_transcripts(sample1_merged_transcripts)

    # print()
    # print()
    _, sample2_merged_transcripts = check_exact_substrings_within_sets(sample2_seqs, args.merge_5prime, args.merge_3prime)
    sample2_merged = merge_transcripts(sample2_merged_transcripts)
    n_s, n_s1, n_s2 = len(final_merged), len(sample1_merged), len(sample2_merged)
    shared = n_s1 + n_s2
    print("Pred_s1,Pred_s2,shared,\%shared")
    if shared > 0:
        print("{0},{1},{2},{3}\n".format(n_s1, n_s2, n_s, round(100*2*n_s/(n_s1 + n_s2), 1)) ) 
    else:
        print("{0},{1},{2},{3}\n".format(n_s1, n_s2, n_s, "-"))
    # n_s1 = len(sample1_seqs)
    # n_s2 = len(sample2_seqs)
    # shared = n_s1 + n_s2
    # print()
    # print()
    # print("Pred_s1,Pred_s2,shared,\%shared")
    # if shared > 0:
    #     print("{0},{1},{2},{3}\n".format(n_s1, n_s2, is_validated, round(100*2*is_validated/(n_s1 + n_s2), 1)))
    # else:
    #     print("{0},{1},{2},{3}\n".format(n_s1, n_s2, is_validated, "-"))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('sample1', type=str, help='Path to consensus fasta file(s)')
    parser.add_argument('sample2', type=str, help='Path to consensus fasta file(s)')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    parser.add_argument('prefix', type=str, help='Prefix file name')

    # parser.add_argument('--names', type=str, nargs=2, required=True, help='Set names')
    # parser.add_argument('--plot_prefix', type=str, required=True, help='Plot file name')
    parser.add_argument('--merge_5prime', type=int, default=100, help='5 prime difference offset')
    parser.add_argument('--merge_3prime', type=int, default=30, help='3 prime difference offset')

    args = parser.parse_args()

    mkdir_p(args.outfolder)

    main(args)





