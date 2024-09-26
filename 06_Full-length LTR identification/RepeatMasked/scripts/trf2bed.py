'''
Descripttion: 
version: 
Author: zpliu
Date: 2023-04-28 11:01:49
LastEditors: zpliu
LastEditTime: 2023-04-28 11:01:49
@param: 
'''
import sys
import pandas as pd
trf = []
header = '#chr start end PeriodSize CopyNumber ConsensusSize PercentMatches PercentIndels Score A C G T Entropy Motif Sequence'.split()
for datf in snakemake.input.dats:
	chrom = None
	sys.stderr.write( "\r" + datf )			
	with open(datf, 'r') as dat:
				for line in dat:
					splitline = line.split()
					if( line.startswith("Sequence:") ):
						chrom = int(line.split()[1].strip())
						#sys.stderr.write(chrom + "\n")
					elif( line.startswith("@") ):
						chrom = splitline[0][1:].strip() # grab everything after the @ in the first word
					else:
						# Catch index errors when line is blank
						try:
							# Check if in header sequence (all non-header lines start with an int: start pos)
							try:
								int(splitline[0])
							except ValueError:
								continue
							trf.append([chrom] + splitline[ 0: (len(header)-1) ] )
						except IndexError:
							pass
trf = pd.DataFrame(trf, columns=header)
print(trf.shape)
		
trf["start"] = trf["start"].astype(int)
trf.sort_values(by=["#chr", "start"], inplace=True)
print("done sorting trf")

trf.to_csv(snakemake.output.bed, sep="\t", index=False)