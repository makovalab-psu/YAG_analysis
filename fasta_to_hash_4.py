#!/usr/bin/env python3
"""
Convert a fasta file to one line per sequence, consisting of the hash code of
the nucleotide sequence, an optional command-line specified identifier, and the
sequence header (including the ">").
"""

from sys  import argv,stdin,stderr,exit
from math import ceil

try:                from hashlib import md5 as md5_new
except ImportError: from md5     import new as md5_new


def usage(s=None):
	message = """
usage: cat fasta_file | fasta_to_hash [<K>/]<N> [<id>] [options]
  --canonical          consider reverse-complemented equivalent sequences to be
                       the same
                       (by default we consider e.g. AAC and GTT to be different)
  [<K>/]<N>            conceptually split the files into <N> groups, and only
                       write the <K>th group; <K> ranges from 1 to <N>;  the
                       numerator can be a comma-separated list, in which case
                       all indicated groups are written;  the numerator can
                       also be a file name, in which case we read a set of
                       groups from that file;  if the numerator is absent, all
                       groups are written
  --complement         output the sequences that are NOT selected by residue
                       (note that this is the mathematical complement of a
                       set, NOT reverse complement of a nucleotide sequence)
  --report=actual      report the actual hash, even though we are filtering
                       with a modulus
  --report=both        report both the modulus and the actual hash
  --report=hashonly    don't report the sequence header
  --report=sequence    also report the nucleotide sequence
  <id>                 optional identifier for output
  --head=<number>      limit the number of input sequences
  --progress=<number>  periodically report how many sequences we've read

  Convert a fasta file to one line per sequence, consisting of the hash code of
  the nucleotide sequence, an optional command-line specified identifier, and
  the sequence header (including the ">").

  Output looks something like this:

	0000223ab87a4735541d30e84 >DS6BPQV01BR6O6
	00040f4759bb6805dd86c146c >DS6BPQV01BHZX9

  or this:

	0000223ab87a4735541d30e84 id >DS6BPQV01BR6O6
	00040f4759bb6805dd86c146c id >DS6BPQV01BHZX9

We use md5 for the hash, but reduce the result to 100 bits. If <N> is
specified, we treat this as a number and reduce it modulo M. Otherwise, we
output it as 25 hex characters."""

	if (s == None): exit (message)
	else:           exit ("%s\n%s" % (s,message))


def main():
	global canonicalize
	global debug

	canonicalize    = False
	N = K           = None
	kFilename       = None
	passNonResidues = False
	reportWhat      = "modulus"
	reportNucs      = False
	id              = ""
	headLimit       = None
	reportProgress  = None
	debug           = []

	for arg in argv[1:]:
		if ("=" in arg):
			argVal = arg.split("=",1)[1]

		if (arg in ["--canonical","--canonicalize","--canon"]):
			canonicalize = True
		elif (arg == "--complement"):
			passNonResidues = True
		elif (arg == "--report=actual"):
			reportWhat = "actual"
		elif (arg == "--report=both"):
			reportWhat = "both"
		elif (arg == "--report=hashonly"):
			reportWhat = "hash only"
		elif (arg == "--report=sequence"):
			reportNucs = True
		elif (arg.startswith("--id=")):
			id = " " + argVal
		elif (arg.startswith("--head=")):
			headLimit = int_with_unit(argVal)
		elif (arg.startswith("--progress=")):
			reportProgress = int_with_unit(argVal)
		elif (arg == "--debug"):
			debug += ["debug"]
		elif (arg.startswith("--debug=")):
			debug += argVal.split(",")
		elif (arg.startswith("--")):
			usage("unrecognized option: %s" % arg)
		elif (N == None):
			if ("/" not in arg):
				N = int_with_unit(arg)
			else:
				(K,N) = arg.split("/",1)
				N = int_with_unit(N)
				if ("," in K):
					K = [int(k) for k in K.split(",")]
					for (ix,k) in enumerate(K):
						assert (k not in K[ix+1:])
						assert (0 < k <= N)
				else:
					try:
						Ki = int(K)
						K = Ki
						assert (0 < K <= N)
						K = [K]
					except ValueError:
						kFilename = K
						K = None
		elif (id == None):
			id = " " + arg
		else:
			usage("unrecognized option: %s" % arg)

	if (passNonResidues) and (K == None) and (kFilename == None):
		usage("you can't use --complement without <K>/<N>")

	# read groups

	if (K == None) and (kFilename != None):
		f = open(kFilename,"rt")
		K = [int(k) for k in f]
		for (ix,k) in enumerate(K):
			assert (k not in K[ix+1:])
			assert (0 < k <= N)
		f.close()

	# process the fasta sequences

	seqNum = 0
	for (seqName,seqHash,seq) in sequences(stdin):
		seqNum += 1
		if (headLimit != None) and (seqNum > headLimit):
			print("limit of %d sequences reached" % headLimit,file=stderr)
			break
		if (reportProgress != None) and (seqNum % reportProgress == 0):
			print("progress: sequence %d" % seqNum,file=stderr)

		if (reportNucs): seqInfo = seqName + " " + seq
		else:            seqInfo = seqName

		if (N == None):
			if (reportWhat == "hash only"): print(             seqHash)
			else:                           print("%s%s %s" % (seqHash,id,seqInfo))
		else:
			filterVal = 1 + md5_to_value(seqHash,N)
			if (K != None):
				if ((filterVal in K) == passNonResidues): continue
			if   (reportWhat == "actual"):    print("%s%s %s"    % (          seqHash,id,seqInfo))
			elif (reportWhat == "both"):      print("%d %s%s %s" % (filterVal,seqHash,id,seqInfo))
			elif (reportWhat == "hash only"): print(                          seqHash)
			else:                             print("%d%s %s"    % (filterVal,        id,seqInfo))


def sequences(f):
	seqName = seqHash = None
	for line in stdin:
		line = line.strip()

		if (line.startswith(">")):
			if (seqName != None):
				seq = seqCanon = "".join(seq).upper()
				if (canonicalize):
					seqRev = reverse_complement(seq)
					seqCanon = min(seq,seqRev)
				seqHash = md5_new()
				seqHash.update(seqCanon.encode("utf-8"))
				yield (seqName,seqHash.hexdigest()[:25],seq)
			(seqName,seq) = (line,[])
		elif (seqName == None):
			assert (False), "first sequence has no header"
		else:
			seq += [line]

	if (seqName != None):
		seq = seqCanon = "".join(seq).upper()
		if (canonicalize):
			seqRev = reverse_complement(seq)
			seqCanon = min(seq,seqRev)
		seqHash = md5_new()
		seqHash.update(seqCanon.encode("utf-8"))
		yield (seqName,seqHash.hexdigest()[:25],seq)


def md5_to_value(s,modulus):
	v = int(s,16) % modulus
	return v


# reverse_complement--

complementMap = str.maketrans("ACGTSWRYMKBDHVNacgtswrymkbdhvn",
                              "TGCASWYRKMVHDBNtgcaswyrkmvhdbn")

def reverse_complement(nukes):
	return nukes[::-1].translate(complementMap)


# int_with_unit--
#	Parse a string as an integer, allowing unit suffixes

def int_with_unit(s):
	if (s.endswith("K")):
		multiplier = 1000
		s = s[:-1]
	elif (s.endswith("M")):
		multiplier = 1000 * 1000
		s = s[:-1]
	elif (s.endswith("G")):
		multiplier = 1000 * 1000 * 1000
		s = s[:-1]
	else:
		multiplier = 1

	try:               return          int(s)   * multiplier
	except ValueError: return int(ceil(float(s) * multiplier))


if __name__ == "__main__": main()
