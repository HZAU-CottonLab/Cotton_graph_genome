import pandas as pd
import numpy as np
import sys

filepath = sys.argv[1]
TEs = []
samples = set()
with open(filepath)as f:
    for index, line in enumerate(f):
        k = set([i.split('-')[0] for i in line.split()[2].split(";")])
        TEs.append(line.split()[0])
        samples = samples | k
samples = list(samples)

samples_num = len(samples)


count = index + 1

# Create an empty dataframe with '0', and the number of rows is the number of LTR matrix.
df = pd.DataFrame(np.zeros([count, samples_num]), columns=samples, )

with open(filepath)as f:
    for index, line in enumerate(f):
        k = set([i.split('-')[0] for i in line.split()[2].split(";")])
        df.loc[[index], k] += 1
df.index = TEs
df.to_csv(f"PAV_{filepath}.csv", sep="\t", )

# Calculate the similarity of LTR between two genomes
s = pd.DataFrame(np.ones([samples_num, samples_num]), index=samples, columns=samples)
for i in range(len(samples)):
    for j in range(i + 1, len(samples)):
        s1 = df.query(f"{samples[i]} == {samples[j]}").shape[0]

        s.iloc[i, j] = s1 / df.shape[0]

        s.iloc[j, i] = s.iloc[i, j]

s.to_csv(f"similarity_{filepath}.txt", sep="\t")

# Generate results for different permutations and combinations of multiple genomes
from itertools import permutations

a = list(permutations(range(samples_num)))
result = []
for i in range(samples_num):
    res = [k[:i + 1] for k in a]
    res = list(set([tuple(sorted(i)) for i in res]))
    result.append(res)
for i in range(len(result)):
    result[i] = [list(i) for i in result[i]]

# Core LTR
num = []
for index, v in enumerate(result):
    for j in v:
        df1 = df.iloc[:, j]
        df1["count"] = df1.sum(axis=1)
        n = df1.query(f"count == {index + 1} ").shape[0]
        num.append([index + 1, n])
num = [str(i) for i in num]
with open(f"{filepath}_Coregene.txt", "w")as f:
    for i in num:
        f.write(i + "\n")
n = []
with open(f"{filepath}_Coregene.txt")as f:
    for i in f:
        i = [int(i) for i in i[1:-2].split(',')]
        n.append(i)
pd.DataFrame(n).to_csv(f"{filepath}_Coregene.txt", sep="\t", index=False, header=None)

# Pan LTR
num = []
for index, v in enumerate(result):
    for j in v:
        df1 = df.iloc[:, j]
        df1["count"] = df1.sum(axis=1)
        n = df1.query(f"count != 0 ").shape[0]
        num.append([index + 1, n])
num = [str(i) for i in num]
with open(f"{filepath}_Pangene.txt", "w")as f:
    for i in num:
        f.write(i + "\n")
n = []
with open(f"{filepath}_Pangene.txt")as f:
    for i in f:
        i = [int(i) for i in i[1:-2].split(',')]
        n.append(i)
pd.DataFrame(n).to_csv(f"{filepath}_Pangene.txt", sep="\t", index=False, header=None)