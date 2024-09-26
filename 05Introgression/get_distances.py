#!/usr/bin/env python
import argparse
import pandas as pd 
import numpy as np


"""
Given a SV vcf file, calculate the maximum Jaccard similarity between each accession and the pimps for non-overlapping
windows along the reference genome coordinates.

Output is a matrix. Each row is an SLL accession, and each column is the non-overlapping window. Each element is 
the maximum Jaccard similarity of SVs in that window between that accession and the pimp accessions. 
"""

parser = argparse.ArgumentParser(description='Get Jaccard similarity between SLL and a comparison group.')
parser.add_argument("-vcf", metavar="<SVs.vcf>", type=str, help="SV vcf file with support vectors. Only one chromosome at a time allowed.")
parser.add_argument("-chr", metavar="<chr_name>", type=str, help="Name of reference chromosome.")
parser.add_argument("-species_file", metavar="<group.txt>", type=str, help="First column is the phylogenetic group (SLC, SP, GAL, CHE, or SLL), second column is the accession")
parser.add_argument("-queryspecies", metavar="<wild>", type=str, default="SP", help="Group to compare to SLL (SLC, SP, GAL, or CHE)")
parser.add_argument("-targetspecies", metavar="<cultivar>", type=str, default="SP", help="Group to compare to SLL (SLC, SP, GAL, or CHE)")
parser.add_argument("-fai", metavar="<reference.fasta.fai>", type=str, help="Fasta index file for the reference genome")
parser.add_argument("-w", metavar="<100000>", type=int, default=1000000, help="Introgression window size.")
parser.add_argument("-m", metavar="5", type=int, default=5, help='minimum number of SVs needed to calculate Jaccard')

args = parser.parse_args()
vcf_file = args.vcf
chr_used = args.chr
fai_file = args.fai
species_file = args.species_file
comp_species = args.queryspecies
target_species=args.targetspecies
window_size = args.w
min_den = args.m

if min_den < 2:
    raise ValueError("-m must be at least 2")

speciesData=pd.read_csv(
    species_file,header=None,sep="\t",index_col=None
)
#* 
targetSpecies=speciesData.loc[speciesData[1]==target_species][0].values
querySpecies=speciesData.loc[speciesData[1]==comp_species][0].values

# Get the chromosome sizes
chr_lens = dict()
with open(fai_file, "r") as f:
    for line in f:
        header, length, x, y, z = line.rstrip().split("\t")
        chr_lens[header] = int(length)

# Iterate through the VCF file and get the distances
# Assumes only one chromosome at a time
acc_vec = None
n_windows = chr_lens[chr_used] // window_size
distances = np.zeros((len(targetSpecies), n_windows))
comp_max_accs = np.zeros((len(targetSpecies), n_windows), dtype=np.int32)
current_window = 0
supp_matrix = []


with open(vcf_file, "r") as f:
    for line in f:
        line = line.rstrip()
        #* 统计行数
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                acc_vec = line.split("\t")[9:]
                assert set(acc_vec) == set(speciesData[0].values)
        elif line.startswith(chr_used):
            fields = line.split("\t")
            tags = fields[7].split(";")
            start = int(fields[1])
            widx = start // window_size

            # Get the support vector
            supp_vec = None
            #* 获取每个样本的基因型数据
            #? 0/0 编码为0，其他情况编码为1
            # for j in tags:
            #     if j.startswith("SUPP_VEC="):
            #         supp_vec = [int(i) for i in j[9:]]
            # if supp_vec is None:
            #     raise ValueError("Missing 'SUPP_VEC' field")
            supp_vec=[ 0 if field.startswith("0/0")  else 1  for field in line.split("\t")[9:]]

            # Build the support matrix for this window or start a new matrix for a new window
            if widx == current_window:
                supp_matrix.append(supp_vec)
            else:
                # We have moved on to the next window. Save the distances for the finished window
                sm = np.asarray(supp_matrix)

                # For now, I will calculate the distances in a 'for' loop. Perhaps vectorize in the future
                # Iterate over the SLLs
                t_distances = [[] for i in range(len(targetSpecies))]
                for i in range(len(targetSpecies)):
                    this_acc = targetSpecies[i]
                    supp_idx = acc_vec.index(this_acc)
                    this_vec = sm[:, supp_idx]

                    # Iterate over the comp species:
                    for comp_acc in querySpecies:
                        comp_supp_idx = acc_vec.index(comp_acc)
                        this_comp_vec = sm[:, comp_supp_idx]
                        if np.count_nonzero(this_comp_vec) >= min_den and np.count_nonzero(this_vec) >= min_den:
                            # Get the Jaccard distance
                            num = np.count_nonzero(np.logical_and(this_vec, this_comp_vec))  # Intersection
                            den = np.count_nonzero(np.logical_or(this_vec, this_comp_vec))  # Union
                            t_distances[i].append(num / den)
                        else:
                            t_distances[i].append(-1)

                # Find which comp sample gave the max
                t_distances_argmax = np.asarray([np.argmax(i) for i in t_distances])
                comp_max_accs[:, current_window] = t_distances_argmax
                # Get the max % shared SVs between a given SLL and each comp species sample
                t_distances = np.asarray([np.max(i) for i in t_distances])

                distances[:, current_window] = t_distances

                # Now that we have calculated the distances for the finished window, start the next one
                if widx == n_windows:
                    #* 最后可能还有几个SVs，但是被打断循环了
                    break
                current_window = widx
                supp_matrix = []
                supp_matrix.append(supp_vec)
        else:
            print(line)
# Write out the comparison species that gave the max
with open("comp_matrix." + comp_species +"-"+chr_used+ ".txt", "w") as f:
    f.write("Sample\t" + "\t".join([str(i*window_size)
            for i in range(n_windows)]) + "\n")
    for i in range(len(targetSpecies)):
        # * 判断当前最相似的donor的similarity是否为-1
        donorList_v1 = []
        for windowIndex, queryIndex in enumerate(comp_max_accs[i, :]):
            if distances[i, windowIndex] == -1:
                # * 添加本身
                donorList_v1.append("NA")
            else:
                donorList_v1.append(querySpecies[queryIndex])
        f.write(targetSpecies[i] + "\t" + "\t".join(donorList_v1) + "\n")
# * 对于每个windows统计每个query 材料分别对于几个target材料

DonorCountArray = np.zeros([querySpecies.shape[0], n_windows],dtype=int)
for windowIndex in range(n_windows):
    maxSimilarityDonor = comp_max_accs[:, windowIndex]
    # * 判断每个target和对应的供体similarity为-1
    for targetIndex, donorIndex in enumerate(maxSimilarityDonor):
        if distances[targetIndex, windowIndex] == -1:
            # * 该target材料没有donor
            pass
        else:
            # * 根据donor index修改对于数据
            DonorCountArray[donorIndex, windowIndex] += 1
DonorCountData = pd.DataFrame(
    DonorCountArray, columns=[str(i*window_size) for i in range(n_windows)],
    index=querySpecies
)
DonorCountData.to_csv("comp_DonorCount.{}-{}.txt".format(comp_species,
                      chr_used), header=True, index=True, sep="\t")



# Print the output matrix
print("Sample\t" + "\t".join( [str(i*window_size) for i in range(n_windows)] ))
for i in range(len(targetSpecies)):
    print(targetSpecies[i] + "\t" + "\t".join( [str(j) for j in list(distances[i, :])] ).replace("-1.0", "NA"))

