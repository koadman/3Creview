#!/usr/bin/env python
# tool to count up proximity ligation rates
# usage: samtools view -q20 | count_plig_rate.py <library name>
from __future__ import division
import sys

if len(sys.argv) < 2:
	sys.exit("Usage: samtools view -q20 | count_plig_rate.py <library name>\n")

close_ins = dict()
plig_ins = dict()
plig_ins_2 = dict()

plig_threshold_1 = 1000
plig_threshold_2 = 50000
cross_species = 0
plasmid_links = dict()
mapped = 0

for line in sys.stdin:
	ll = line.split()
	ins = int(ll[8])
	if(ll[6] != "=" and ll[6] != "*"):
		if(ll[2] == "gi|685631213|gb|CP009362.1|"):
			plasmid_links[ll[6]] = plasmid_links.get(ll[6],0) + 1

		cross_species += 1
		continue 
	if(ins < 1):
		continue	
	mapped += 1
	if(ins < plig_threshold_1):
		close_ins[ll[2]] = close_ins.get(ll[2],0) + 1
		continue
	plig_ins[ll[2]] = plig_ins.get(ll[2],0) + 1	
	if(ins < plig_threshold_2):
		plig_ins_2[ll[2]] = plig_ins_2.get(ll[2],0) + 1	


for org,count in close_ins.items():
	thisorg = plig_ins.get(org,0) + close_ins.get(org,1)
	ratio = plig_ins.get(org,0) / thisorg
	ratio_2 = plig_ins_2.get(org,0) / thisorg

	print sys.argv[1] + "\t" + org + "\t" + str(ratio) + "\t" + str(ratio_2) + "\t" + str(thisorg / mapped)

sys.stderr.write("x-species\t" + str(cross_species) + "/" + str(mapped) + " = " + str(cross_species / mapped) + "\n")

for org,count in plasmid_links.items():
	sys.stderr.write("plasmid to " + org + "\t" + str(count) + "\n")
