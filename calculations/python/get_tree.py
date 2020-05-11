from Bio import Phylo
import subprocess, networkx, re
from functools import reduce
from calculations.python.paths import *


#Phylo incorrectly parses the gisaid_china.MCC.trees file, so we modify it, assume that EPI_ISL_\d+ is an identifier
subprocess.check_output(f'sed -E "s/h\S+EPI_ISL_/EPI_ISL_/" {data_path(TREE_UNFILTERED)} | sed -E "s/\|\S+[^,]//" > {data_path(TREE_FILTERED)}', shell=True)
trees = Phylo.parse(data_path(TREE_FILTERED), 'nexus')
tree = trees.__next__()
net = Phylo.to_networkx(tree).to_undirected()

#descendants_at_level[i] keeps a dict {node: nodes that can be reached with exactly i steps}
descendants_at_level = [{node: {node} for node in net.nodes}]

LEVEL = 2

for i in range(LEVEL):
    prev_level = descendants_at_level[i]
    descendants_at_level.append({node: reduce(set.union, (set(net.neighbors(n)) for n in prev_level[node])) for node in net.nodes})

#dict actually used to generate list of genomes to compare
d = {node: reduce(set.union, (descendants_at_level[i][node] for i in range(LEVEL + 1))) for node in net.nodes}

pairs = {frozenset({node.name, descendant.name}) for node in net.nodes for descendant in d[node] if node.name is not None and descendant.name is not None and node.name != descendant.name}

pairs = [tuple(p) for p in pairs]

names = subprocess.check_output(f"cat {data_path(FILTERED_GENOMES)} | grep EPI_ISL_", shell=True).split()
names = [name.decode("UTF-8") for name in names]
names = {re.search("EPI_ISL_\\d+", name).group(): name for name in names}

new_pairs = [(a, b) for a, b in pairs if a in names and b in names]

print(f"{len(pairs) - len(new_pairs)} pairs removed, {len(new_pairs)} left due to missing names.")

file = open(data_path(EXTRA_COMPARISONS_REQUESTS), "w")
file.write(str(len(new_pairs)) + "\n")
for a, b in new_pairs:
    file.write(names[a])
    file.write("\n")
    file.write(names[b])
    file.write("\n")
file.close()
