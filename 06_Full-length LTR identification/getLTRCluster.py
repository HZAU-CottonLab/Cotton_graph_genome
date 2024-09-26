#!/usr/bin/env python3
import sys
import re

file = sys.argv[1]
hash = {}
id = ""

with open(file, 'r') as f:
    for line in f:
        line = line.strip()
        if re.match(r'^>Cluster\s(.*)', line):
            id = "Cluster" + re.match(r'^>Cluster\s(.*)', line).group(1)
        else:
            regions = line.split('\t')
            match = re.search(r'>(.*)\.\.\. at', regions[1])
            if match:
                name = match.group(1)
                if id not in hash:
                    hash[id] = []
                hash[id].append(name)

for t, TEs in hash.items():
    num = len(TEs)
    print(t, num, ';'.join(TEs), sep='\t')