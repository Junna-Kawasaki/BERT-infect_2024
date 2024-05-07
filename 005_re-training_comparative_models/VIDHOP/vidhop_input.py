#!/usr/bin/env python
import sys, re
argvs = sys.argv

import os
import collections
import numpy as np
import pandas as pd

## argvs
fasta_f = argvs[1]
metadata_f = argvs[2]
output_seq_f = argvs[3]
output_host_f = argvs[4]

## read fasta file
fasta_l = []
for line in open(fasta_f):
	line = line.strip()
	if ">" in line:
		query = line.replace(">", "").split("_sliding")[0]
		fasta_l.append([query, ""])
	else:
		seq = line
		fasta_l[-1][1] += seq

df_fasta = pd.DataFrame(fasta_l, columns=["refseq.id", "sequences"])

## read meatadata
df_metadata = pd.read_csv(metadata_f)

## merge df_fasta and df_metadata
df_fasta_host = df_fasta.merge(df_metadata.loc[:, ["refseq.id", "host_label"]],
					on="refseq.id", how="left")

## output
df_fasta_host.to_csv(output_seq_f, sep="\n", columns=["sequences"], 
						index=False, header=False)

df_fasta_host.to_csv(output_host_f, sep="\n", columns=["host_label"],
						index=False, header=False)