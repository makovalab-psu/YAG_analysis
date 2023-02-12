from __future__ import print_function
import os,sys
import argparse
import re
import math
import errno

import re
from pathlib import Path
from collections import defaultdict


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


def main(args):


    directory_orfs = os.fsencode(args.orf_folder)
    all_orf_dict = {}
    for file in os.listdir(directory_orfs):
        filename = os.fsdecode(file)
        info = filename.split('_50aa')[0]
        gene_fam = re.split('\d+_', info)[1]
        exp_id = info.split('_')[0]
        # print(exp_id, gene_fam, filename)

        if exp_id not in all_orf_dict:
            all_orf_dict[exp_id] = {}
        if gene_fam not in all_orf_dict[exp_id]:
            all_orf_dict[exp_id][gene_fam] = set()

        orfs = {'_'.join(acc.split('_')[:-1]) : seq for acc, (seq,_) in readfq(open(os.path.join(args.orf_folder, filename), 'r'))}
        for acc in orfs:
            all_orf_dict[exp_id][gene_fam].add(acc)
        # print({acc for acc in orfs})

    # print(all_orf_dict)
    print()
    print()
    all_transcript_dict = {}

    pathlist = Path(args.shared_transcript_folder).rglob('*.fa')
    for path in pathlist:
        path_in_str = str(path)
        # print(path_in_str.split('/'))
        _, _, _, _, _, _, exp_id, gene_fam, _ = path_in_str.split('/')
        # print(exp_id, gene_fam, path_in_str)
        if exp_id not in all_transcript_dict:
            all_transcript_dict[exp_id] = {}
        if gene_fam not in all_transcript_dict[exp_id]:
            all_transcript_dict[exp_id][gene_fam] = set()
        orfs = {acc: seq for acc, (seq,_) in readfq(open(path_in_str, 'r'))}
        for acc in orfs:
            all_transcript_dict[exp_id][gene_fam].add(acc)

    # print(all_transcript_dict)
    no_orf = defaultdict(lambda: defaultdict(int))
    no_orf_pry = defaultdict(set)

    for exp_id in all_transcript_dict:
        for gene_fam in all_transcript_dict[exp_id]:
            # print(exp_id, gene_fam)
            for tr_acc in all_transcript_dict[exp_id][gene_fam]:
                if tr_acc not in all_orf_dict[exp_id][gene_fam]:
                    # print('NO ORF', tr_acc, exp_id, gene_fam)
                    no_orf[exp_id][gene_fam] += 1
                    if gene_fam == 'PRY':
                        no_orf_pry[exp_id].add(tr_acc)
                # else:
                #     print('HAD ORF', tr_acc, exp_id, gene_fam)

    # print()
    # print()
    for exp_id in no_orf:
        for gene_fam in no_orf[exp_id]:
            print(exp_id,gene_fam,no_orf[exp_id][gene_fam])

    for exp_id in no_orf_pry:
        for tr in no_orf_pry[exp_id]:
            print('{0},{1}'.format(exp_id,tr))
        print()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('orf_folder', type=str, help='Path to files for predicted ORFS')
    parser.add_argument('shared_transcript_folder', type=str, help='Path to files for shared transcripts')

    args = parser.parse_args()

    main(args)
