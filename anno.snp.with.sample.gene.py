#!/usr/bin/env python3
#coding=utf-8
import sys

if len(sys.argv) != 7:
	print ('python %s <clean.snp.pileup2> <sample name> <038.pileup2.snp> <038.glimmer.gff> <vfdb.anno.xls> <prefix>'  % sys.argv[0])
	sys.exit()

from collections import defaultdict
import numpy as np
import re
import pandas as pd

all_snp, sample_name, sample_snp, sample_gff, sample_anno, prefix = sys.argv[1:]
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


snp = AutoVivification()
snp_index = AutoVivification()
anno = AutoVivification()

all_snp_h = open(all_snp, 'r')
header = all_snp_h.readline().strip().split('\t')
header.insert(0, '')
print ('header is:' + ' '.join(header) + ';first e is:' + header[0])

## read vfdb anno

sample_anno_h = open(sample_anno, 'r')
sample_anno_header =sample_anno_h.readline().strip().split('\t') 
sample_anno_header.insert(0, '')
for line in sample_anno_h:
    line = line.strip()
    arr = line.split('\t')
    anno['vfdb'][sample_name][arr[0]]['Description'] = arr[6]
    anno['vfdb'][sample_name][arr[0]]['VF_id'] = arr[3]
    anno['vfdb'][sample_name][arr[0]]['GI_id'] = arr[4]
    anno['vfdb'][sample_name][arr[0]]['Type'] = arr[5]

sample_anno_h.close()


## read gff file
gff = AutoVivification()
sample_gff_h = open(sample_gff, 'r')
for line in sample_gff_h:
    if re.match("#", line) : continue
    line = line.strip()
    arr = line.split('\t')
    if arr[2] == 'gene' :
        re.compile(r'ID=[^;]+').finditer(arr[-1])
        geneid = re.findall('ID=([^;]+)', arr[-1])[0]
        gff[arr[0]][geneid]['start'] = int(arr[3])
        gff[arr[0]][geneid]['end'] = int(arr[4])
        gff[arr[0]][geneid]['strand'] = arr[6]
sample_gff_h.close()

def check_sample_gene(scf, pos):
    for gene in gff[scf]:
        if pos >= gff[scf][gene]['start'] and pos <= gff[scf][gene]['end']:
            return gene

##read all snp
for line in all_snp_h:
	line = line.strip()
	arr = line.split('\t')
	for index in np.arange(2, len(header)):
		if arr[index] != arr[1]:
                        #snp[header[index]][arr[0]]["sample"] = arr[index]
                        snp[header[index]][arr[0]]['sample'] = arr[index]
                        snp[header[index]][arr[0]]['ref'] = arr[1]

all_snp_h.close()

## read single sample ref2scf
sample_snp_h = open(sample_snp, 'r')
for line in sample_snp_h:
	line = line.strip()
	arr = line.split('\t')
	snp_index_tmp = arr[4] + '_' + arr[0]
	if snp_index_tmp in snp[sample_name]:
                scaffold = arr[5]
                scaffold_pos = arr[2]
                sample_gene = check_sample_gene(scaffold, int(scaffold_pos))
                snp[sample_name][snp_index_tmp]['pos'] = scaffold_pos
                snp[sample_name][snp_index_tmp]['scaffold'] = scaffold
                snp[sample_name][snp_index_tmp]['pos_strand'] = arr[-1]
                #Identity   E_value VF_id   GI_id   Type    
                if sample_gene in anno['vfdb'][sample_name]:
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_Description'] = anno['vfdb'][sample_name][sample_gene]['Description']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_VF_id'] = anno['vfdb'][sample_name][sample_gene]['VF_id']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_GI_id'] = anno['vfdb'][sample_name][sample_gene]['GI_id']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_Type'] = anno['vfdb'][sample_name][sample_gene]['Type']
                else:
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_Description'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_VF_id'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_GI_id'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb_Type'] = "null"


                if sample_gene:
                    snp[sample_name][snp_index_tmp]['sample_gene_id'] = sample_gene
                    snp[sample_name][snp_index_tmp]['sample_gene_id_start'] = gff[scaffold][sample_gene]['start']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_end'] = gff[scaffold][sample_gene]['end']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_strand'] = gff[scaffold][sample_gene]['strand']
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb'] = "tmp"


                else:
                    sample_gene = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_start'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_end'] = "null"
                    snp[sample_name][snp_index_tmp]['sample_gene_id_vfdb'] = "null"


                print ('gene:' + scaffold + ':' + scaffold_pos + '::' + sample_gene + ':::' + line)
                
	

pd.DataFrame(snp[sample_name], index=('pos_strand', 'pos', 'ref',  'sample', 'sample_gene_id', 'scaffold', 'sample_gene_id_start', 'sample_gene_id_end', 'sample_gene_id_strand', 'sample_gene_id_vfdb_Description', 'sample_gene_id_vfdb_VF_id', 'sample_gene_id_vfdb_GI_id', 'sample_gene_id_vfdb_Type')).T.to_csv(prefix + ".out", sep='\t')

