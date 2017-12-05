#!/usr/bin/env python3
#coding=utf-8
import sys

if len(sys.argv) != 3:
	print ('python %s <*gff> <out> '  % sys.argv[0])
	sys.exit()

from collections import defaultdict
import re

gff, out = sys.argv[1:]
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


#gff = AutoVivification()
out_content = ''
out_exon = []
gff_h = open(gff, "r")
for line in gff_h:
	line = line.strip()
	if re.match("#", line):
		out_content = out_content + line + "\n"
		continue
	arr = line.split('\t')
	if arr[2] == "mRNA":
		if len(out_exon) != 0:
			out_content = out_content + "\n".join(out_exon) + "\n" + line + "\n"
			out_exon = []
		else:
			out_content = out_content + line + "\n"

	elif arr[2] == "CDS":
		if arr[6] == "+":
			out_exon.append(line)
		elif arr[6] == "-":
			out_exon.insert(0, line)
		else:
			print ("error")
			sys.exit()
open(out, "w").write(out_content)
