#!/usr/bin/env python

import sys, re
argvs = sys.argv

error1 = "no argvs"
try:
	for i in range(1,10):
		if argvs[i] == "-f":
			fasta_f = argvs[i+1]
		elif argvs[i] == "-n":
			criteria = argvs[i+1]
			gap_count = "n"
		elif argvs[i] == "-p":
			criteria = argvs[i+1]
			gap_count = "p"
except:
	try:
		if len(fasta_f) != 0 and len(criteria) != 0:
			del error1
	except:
		print (error1)

fasta_d = {}
for line in open(fasta_f):
	line = line.strip()
	if ">" in line:
		name = line
		#print (name)
		seq = ""
	else:
		seq = seq + "".join(map(str,line))
		fasta_d[name] = seq

out_d = {}
for name,seq in fasta_d.items():
	seq = [s for s in seq]
	gap_n = seq.count("-")
	gap_p = gap_n / len(seq) * 100
	#print (name,gap_n)
	if gap_count == "n":
		if int(gap_n) < int(criteria):
			out_d[name] = seq
			
	elif gap_count == "p":
		if int(gap_p) < int(criteria):
			out_d[name] = seq

print (len(fasta_d))
print (len(out_d))
print ("\n".join(out_d.keys()))

out = ".".join(map(str,re.split("\.",fasta_f)[:-1])) + "_rmgap.fasta"
f = open(out, "w")
for name,seq in out_d.items():
	f.write(name + "\n")
	f.write("".join(map(str,seq)) + "\n")
f.close()