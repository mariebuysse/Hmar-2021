#!/usr/bin/env python3
import sys
import re
# Use as: python3 split_cogtab.py orthogroups.tsv <GENOME>
genome = sys.argv[2]
orthofile = sys.argv[1]
f = open(orthofile)
for line in f:
	linea=line.strip().split("\t")
	for l in linea[1:]:
		if l == genome:
			print (linea[0])
