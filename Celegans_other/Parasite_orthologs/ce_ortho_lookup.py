#!/bin/python
#python ce_ortho_lookup.py ceQTLvars.csv c_elegans.PRJNA13758.WS255.orthologs.txt > ceQTLvars-final.csv
import collections
from collections import Counter
from sys import argv 
script, filename1, filename2 = argv 
import re
import math
import numpy as np


orthotxt = []
b = open(filename2) #open c_elegans.PRJNA13758.WS255.orthologs.txt 
for line in b:
	line = line.strip()
	orthotxt.append(line)
orthotxt = '\t'.join(orthotxt)
orthotxt = orthotxt.split("=")


a = open(filename1) #open ceQTLvars.csv
for line in a:
	currentline = line.strip()
	currentline = currentline.split(",")
	currentID = currentline[3]
	currentline.append("\"NA\"")
	currentline.append("\"NA\"")
	#print ",".join(currentline)
	currentID = re.sub('\"','', currentID)
	for ortholine in orthotxt:
		m = re.search(currentID,ortholine)
		if m:
			brugia_hits = re.findall('Brugia malayi\t(WBGene\d+)',ortholine)
			if len(brugia_hits) > 0:
				brugia_entry = "; ".join(brugia_hits)
			else:
				brugia_entry = "NA"
			sratti_hits = re.findall('Strongyloides ratti\t(WBGene\d+)',ortholine)
			if len(sratti_hits) > 0:
				sratti_entry = "; ".join(sratti_hits)
			else:
				sratti_entry = "NA"
			currentline[6] = "\""+brugia_entry+"\""
			currentline[7] = "\""+sratti_entry+"\""		
			#currentline.append("\""+brugia_entry+"\"")
			#currentline.append("\""+sratti_entry+"\"")
	print ",".join(currentline)