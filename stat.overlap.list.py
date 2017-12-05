#!/usr/bin/env python3
import sys
if len(sys.argv)!= 3 :
	print (f'python {sys.argv[0]} <new,list> <old.list>\n')
	sys.exit()

new, old = sys.argv[1:]


def fetch_id(file):
	file_hash = {}
	file_h = open(file, "r")
	for line in file_h:
		arr = line.split("\t")
		file_hash[arr[0]] = ""
	return file_hash


new_ids = fetch_id(new)
old_ids = fetch_id(old)

print (f'{old} has,but {new} not:\n')
for old_id in old_ids:
	if old_id not in new_ids:
		print (f'\t{old_id}')

print (f'\n{new} has,but {old} not:\n')
for new_id in new_ids:
	if new_id not in old_ids:
		print (f'\t{new_id}')


