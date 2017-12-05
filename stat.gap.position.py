#!/usr/bin/env /ifshk7/BC_PS/sikaiwei/python/hk_new/Python-2.7.13/bin/python

import sys
if len(sys.argv) < 2:
    print 'python %s <genome.fa> <outfile>' % sys.argv[0]
    sys.exit()
from Bio import SeqIO
import re

genome = sys.argv[1]
outfile = sys.argv[2]
outfile_handle = open(outfile,'w+')
for seq_record in SeqIO.parse(genome,"fasta"):
    id = re.sub('\s+.*$','',seq_record.id)
    fasta = str(seq_record.seq)
    search = re.compile(r'N+').finditer(fasta) ## find all ,search or match only once
    for e in search:
        outfile_handle.write("%s\t%s\t%s\n" % (id, e.start(), e.end()))

outfile_handle.close()

