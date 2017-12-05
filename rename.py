#!/usr/bin/env python

import sys
import re
from collections import defaultdict
from Bio import SeqIO

m8_dict = defaultdict(dict)

m8_file = "a-b.blastn.m8"
new2out = open("new2did","w")
m8_handle = open(m8_file,"r")
for line in m8_handle:
	arr = line.split("\t")
	if arr[2] < 90 or arr[3] < 400:	continue
	m8_dict[arr[0]][arr[1]] = m8_dict[arr[0]].setdefault(arr[1],0) + abs(int(arr[6]) - int(arr[7]))




genome = "adjust.strand.rename.fa"
seq = {}
for record in SeqIO.parse(genome,"fasta"): seq[record.id] = str(record.seq)


id2id = open("GCF_000146045.2_R64_genomic.fna.id.2.id","r")
id2id_d = {}
out = {}
for line in id2id:
	arr = re.split('\s+',line.strip())
	id2id_d[arr[0]] = " ".join(arr[1:])
for ctg in m8_dict:
	ref_name = max(m8_dict[ctg], key=m8_dict[ctg].get)
	new2out.write(ctg + "\t" + ref_name + "\t" + id2id_d[ref_name] + "\n")
	ctg = re.sub('\|.*','',ctg)
	out[id2id_d[ref_name]] = seq[ctg].strip()

new2out.close()
for ctg in sorted(out):
	ctg2 = re.sub("\s+","_",ctg)
	print ">%s\n%s" % (ctg2,out[ctg])



