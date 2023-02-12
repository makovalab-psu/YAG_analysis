
from __future__ import print_function
import os,sys
import argparse
import re


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

def main(params):

    outfile = open(params.outfile, "w")
    p = r'ATG+'
    for acc, (seq, qual) in  readfq(open(params.fastq_file, 'r')):
        m = re.match(p, seq)
        seq_mod = seq
        if m:
            mm = m.group()
            seq_mod = seq[len(mm):]
            # print(mm, len(mm), seq_mod[:20], seq_mod[len(mm):20] )

        outfile.write("@{0}\n{1}\n{2}\n{3}\n".format(acc, seq_mod, "+", qual))



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Evaluate pacbio IsoSeq transcripts.")
    parser.add_argument('fastq_file', type=str, help='Path to fastq file')
    parser.add_argument('outfile', type=str, help='Output path of results')

    args = parser.parse_args()
    main(args)