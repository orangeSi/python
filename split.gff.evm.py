#!/usr/bin/env /ifshk7/BC_PS/sikaiwei/python/Python-2.7.9/bin/python
# -*- coding: UTF-8 -*-
import sys
if len(sys.argv) < 7:
	print "python ",sys.argv[0]," <outdir> <prefix> <genome.fasta> <weights.txt> <genemarkes.gff> <august.gff> <...gff>, at least two gff file\n"
	sys.exit()

import re
import os
from Bio import SeqIO
#from  collections  import defaultdict
#genes = defaultdict(defaultdict)
genes = {}
outdir = os.path.abspath(sys.argv[1])
prefix = sys.argv[2]
genome = os.path.abspath(sys.argv[3])
weight = os.path.abspath(sys.argv[4])

for gff in sys.argv[5:]:
	print '#start read ',gff
	file = open(gff,"r")
	for line in file.readlines():
		if re.search(r'^#',line): continue # if in one line
		arr = line.split("\t")
		genes.setdefault(arr[0],{}) #defined defalut duowei dict
		genes[arr[0]].setdefault(arr[1],"") #defined defalut duowei dict

		genes[arr[0]][arr[1]] += line
		
		line = line.strip()
		#print 'this line is:',line
#	line 13 end
	print '#end read ',gff,'\n'
	
#line 10 end

if outdir != "" and not os.path.exists(outdir): os.makedirs(outdir) 
process = outdir+"/process/"
if not os.path.exists(process): os.makedirs(process) 
shell = outdir+"/shell/"
if not os.path.exists(shell): os.makedirs(shell) 


## read genome file
genomes = {}
for seq_record in SeqIO.parse(genome,"fasta"):
	id = re.sub('\s.*','',seq_record.id)
	genomes[str(id)] = str(seq_record.seq)


## split gff
gff = {}
for scf in genes:
	out_gff = process + prefix + '.' + scf + '.gff'
	gff[scf] = out_gff
	the_scaf_file = open(out_gff,"w+")
	the_scf = ''
	print '#open file ' + out_gff
	for method in genes[scf]: the_scf += genes[scf][method]
	the_scaf_file.write(the_scf)
	the_scaf_file.close()
	print '#close file ' + out_gff
	if scf in genomes:
		out_scf = process + prefix + '.' + scf + '.fasta'
		seq_record_file = open(out_scf,"w+")
		seq_record_file.write('>'+ scf + '\n' + genomes[scf])
		seq_record_file.close()
		print 'write done:' + out_scf
		out_shell_file = shell + prefix + '.' + scf + '.evm.sh'
		out_shell = '''. /ifshk7/BC_PS/sikaiwei/software/EVidenceModeler-1.1.1/env.sh
		## run EVM
		out=%s
		evidence_modeler.pl --genome  %s --weights %s --gene_predictions $out > $out.evm.out && echo Still_waters_run_deep > %s.sign
		for j in $(cat %s|awk '{print $1}'|sort -u)
		do
			EVM_to_GFF3.pl $out.evm.out $j > $out.evm.out.gff3
		done
		'''
		out_shell_file_handle = open(out_shell_file,"w+")
		out_shell_file_handle.write(out_shell % (out_gff,out_scf,weight,out_shell_file,out_gff))
		out_shell_file_handle.close()
		#print out_shell % (out_gff,out_scf,weight,out_gff,out_gff)

	else:
		print 'die:' + id + ' is not in ' + gff
		sys.exit()



