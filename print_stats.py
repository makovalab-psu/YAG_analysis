
from __future__ import print_function
import os,sys
import argparse
from collections import defaultdict
import errno


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



def mkdir_p(path):
    # print("creating", path)
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def read_tsv(f):
    tsv_info = defaultdict(list)
    for line in open(f, 'r'):
        read_acc, ref_acc, read_len, ref_len = line.split()
        tsv_info[ref_acc].append((read_acc, read_len, ref_len))

    return tsv_info


def main(args):
    transcript_name = os.path.basename(args.final_fa_file).split(".fa")[0]
    # print(transcript_name)
    final_fa_seqs = [ (acc, seq.upper()) for acc, (seq, _) in  readfq(open(args.final_fa_file, 'r'))]
    repl1_seqs = { acc : seq.upper() for acc, (seq, _) in  readfq(open(args.reads_repl1, 'r'))}
    repl2_seqs = { acc : seq.upper() for acc, (seq, _) in  readfq(open(args.reads_repl2, 'r'))}
    tsv_info1 = read_tsv(args.tsv_repl1)
    tsv_info2 = read_tsv(args.tsv_repl2)
    # print(list(tsv_info1.keys()))
    # csv_info_file = open(final_fa_file + "_INFO.csv", "w")
    f = open(args.outfolder + "/" + transcript_name + "_reads.fa", "w")
    transcript_length = -1
    read_lengths = []
    supporting_reads_fa = mkdir_p(args.outfolder)
    for i, (acc, consensus_seq) in enumerate(final_fa_seqs):
        if i == 0:
            transcript_length = len(consensus_seq)
        if "REPLICATE1" in acc and "REPLICATE2" in acc:
            acc_repl1 = "_".join(acc.split("_")[1:5])
            # print("HERE", acc_repl1)
            # print(acc.split("_"))
            # acc_repl1 = acc.split("REPLCATE2")[0].split("CATE1_")[1]
            acc_repl2 = acc.split("REPLICATE2_")[1]
            acc_repl2 = "_".join(acc_repl2.split("_")[:4])
            for (read_acc, read_len, ref_len) in tsv_info1[acc_repl1]:
                seq = repl1_seqs[read_acc]
                f.write(">{0}\n{1}\n".format(read_acc, seq))
                read_lengths.append(len(seq))
            # print(nr_reads)
            for (read_acc, read_len, ref_len) in tsv_info2[acc_repl2]:
                seq = repl2_seqs[read_acc]
                f.write(">{0}\n{1}\n".format(read_acc, seq))
                read_lengths.append(len(seq))
            # print(nr_reads)
        elif "REPLICATE2_" in acc:
            acc_repl2 = "_".join(acc.split("_")[1:5])
            for (read_acc, read_len, ref_len) in tsv_info2[acc_repl2]:
                seq = repl2_seqs[read_acc]
                f.write(">{0}\n{1}\n".format(read_acc, seq))
                read_lengths.append(len(seq))

        elif "REPLICATE1_" in acc:
            acc_repl1 = "_".join(acc.split("_")[1:5])
            for (read_acc, read_len, ref_len) in tsv_info1[acc_repl1]:
                seq = repl1_seqs[read_acc]
                f.write(">{0}\n{1}\n".format(read_acc, seq))
                read_lengths.append(len(seq))

    r = sorted(read_lengths)
    nr_reads = len(r)
    min_read_len = min(r)
    if len(r) % 2 == 1:
        m = len(r)//2
        median_read_len = r[m]
    else:
        m1 = len(r)//2
        m2 = m1 - 1
        v1 = r[m1]
        v2 = r[m2]
        median_read_len = (v1 + v2)/2
    max_read_len = max(r)
    # for the final INFO CSV file
    # transcript name, transcript length, nr_reads, min_read_len, median_read_len, max_read_len
    print("{0},{1},{2},{3},{4},{5}".format(transcript_name, transcript_length, nr_reads, min_read_len, median_read_len, max_read_len) )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('final_fa_file', type=str, help='Path to consensus fasta file')
    parser.add_argument('reads_repl1', type=str, help='Path to reads replicate 1')
    parser.add_argument('reads_repl2', type=str, help='Path to reads replicate 2')
    parser.add_argument('tsv_repl1', type=str, help='Path to IsoCon tsv file replicate 1')
    parser.add_argument('tsv_repl2', type=str, help='Path to IsoCon tsv file replicate 2')
    parser.add_argument('outfolder', type=str, help='Output path of results')
    # parser.add_argument('prefix', type=str, help='Prefix file name')
    args = parser.parse_args()

    mkdir_p(args.outfolder)

    main(args)





