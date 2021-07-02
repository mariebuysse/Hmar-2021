#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 10:32:57 2019
"""
from Bio import SeqIO
import sys
with open(sys.argv[2]) as c:
	lista_nomi=c.read().splitlines()
out = (str(sys.argv[1])+".subset.fasta")
if len(sys.argv)==4:
		if sys.argv[3]=="-r":
			print("Excluding the list")
			records = (r for r in SeqIO.parse(sys.argv[1], "fasta") if r.id not in lista_nomi)
		else:
			print("Wrong sysargv")
else:
	records = (r for r in SeqIO.parse(sys.argv[1], "fasta") if r.id in lista_nomi)
SeqIO.write(records, out, "fasta")
