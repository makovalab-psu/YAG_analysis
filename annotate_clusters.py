
from __future__ import print_function
import os,sys
import argparse
import re
import math
import edlib
import errno

def reverse_complement(string):
    #rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N', 'X':'X'}
    # Modified for Abyss output
    rev_nuc = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'a':'t', 'c':'g', 'g':'c', 't':'a', 'N':'N', 'X':'X', 'n':'n', 'Y':'R', 'R':'Y', 'K':'M', 'M':'K', 'S':'S', 'W':'W', 'B':'V', 'V':'B', 'H':'D', 'D':'H', 'y':'r', 'r':'y', 'k':'m', 'm':'k', 's':'s', 'w':'w', 'b':'v', 'v':'b', 'h':'d', 'd':'h'}

    rev_comp = ''.join([rev_nuc[nucl] for nucl in reversed(string)])
    return(rev_comp)

def edlib_traceback(x, y, additional_equalities, mode="HW", task="path", k=1):
    result = edlib.align(x, y, mode=mode, task=task, k=k, additionalEqualities= additional_equalities)
    ed = result["editDistance"]
    if ed == -1:
        ed = len(x)

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
    probes = {seq.upper() : acc for acc, (seq, _) in  readfq(open(args.probes, 'r'))}
    additional_equalities = [("R", "A"), ("R", "G"), 
                            ("Y", "C"), ("Y", "T"),
                            ("M", "A"), ("M", "C"),
                            ("K", "G"), ("K", "T"),
                            ("S", "C"), ("S", "G"),
                            ("W", "A"), ("W", "T"),
                            ("H", "A"), ("H", "C"), ("H", "T"),
                            ("B", "C"), ("B", "G"), ("B", "T"),
                            ("V", "A"), ("V", "C"), ("V", "G"),
                            ("D", "A"), ("D", "G"), ("D", "T"),
                            ("N", "A"), ("N", "C"), ("B", "G"), ("B", "T")]
    # additional_equalities = []
    fam_counts = {}
    for read_acc, (read_seq, _) in  readfq(open(args.reads, 'r')):
        best_ed = 100
        best_acc = "NA"
        for probe_seq, gene_fam_acc in probes.items():
            # probe is query, read is target
            (ed, locations, cigar) = edlib_traceback(probe_seq, read_seq, additional_equalities, task="path", k=best_ed)
            # print(gene_fam_acc, "ed:", best_ed, "to:", best_acc)
            if ed < best_ed:
                best_ed = ed
                best_acc = gene_fam_acc

            read_seq_rc = reverse_complement(read_seq)
            (ed, locations, cigar)  =  edlib_traceback(probe_seq, read_seq_rc, additional_equalities, task="path", k=best_ed)
            # print(gene_fam_acc, "RC ed:", best_ed, "to:", best_acc)
            if ed < best_ed:
                best_ed = ed
                best_acc = gene_fam_acc

        print("Best ed:", best_ed, "to:", best_acc, cigar, read_acc)
        if best_acc in fam_counts:
            fam_counts[best_acc] += 1/max(1, best_ed)
        else:
            fam_counts[best_acc] = 1/max(1, best_ed)
    print(fam_counts)
    if fam_counts:
        file_name_gen_fam = sorted([(count, acc) for acc, count in fam_counts.items()], key= lambda x: x[0], reverse=True)[0][1]
    else:
        file_name_gen_fam = "unannotated"
    print(file_name_gen_fam)
    p = args.prefix + "_" + file_name_gen_fam + ".fa"
    outfile = open(os.path.join(args.outfolder, p), "a")
    for read_acc, (read_seq, _) in  readfq(open(args.reads, 'r')):
        outfile.write(">{0}\n{1}\n".format(read_acc, read_seq))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Annotates a cluster file with the best matching header of probe sequences.")
    parser.add_argument('reads', type=str, help='reads)')
    parser.add_argument('probes', type=str, help='probes')
    parser.add_argument('outfolder', type=str, help='Output folder')
    parser.add_argument('prefix', type=str, help='Outfile prefix')

    args = parser.parse_args()

    mkdir_p(args.outfolder)

    main(args)





