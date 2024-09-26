import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load the data
df = pd.read_table("AtDt_gene_presence_absence.txt")
data = df.iloc[:, 1:].astype(np.int8)

core_num = []
pan_num = []

# Define the number of iterations
num_iterations = 10000

for index in range(1, 36):  # Number of samples + 1
    core = []
    pan = []
    for _ in range(num_iterations):
        sampled_data = data.sample(index, axis=1, replace=False).sum(axis=1)
        core.append((sampled_data == index).sum())
        pan.append((sampled_data != 0).sum())
    core_num.append(core)
    pan_num.append(pan)
    print(f"Processed {index} samples")

# Save core genome data to file
with open("core_num.txt", "w") as f:
    for index, values in enumerate(core_num):
        f.write(f"{index + 1}\t" + "\t".join(map(str, values)) + "\n")

# Save pan genome data to file
with open("pan_num.txt", "w") as f:
    for index, values in enumerate(pan_num):
        f.write(f"{index + 1}\t" + "\t".join(map(str, values)) + "\n")