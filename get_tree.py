from Bio import Phylo

import networkx, pylab, subprocess, re

tree = Phylo.read('tree2.nwk', 'newick')

net = Phylo.to_networkx(tree)

list = networkx.generate_adjlist(net)
#print(net.clade_relations)
#print(list)
counter = 0

import sys
sys.setrecursionlimit(150000)
def FindNeighboors(list, node, counter):
    for x in net.adj[node]:
        if x.name:
            list.append(x.name)
        else:
            if counter<100:
                print(counter)
                counter += 1
                FindNeighboors(list,x,counter)
            else:break
    return list, node,counter #that wont work

def FindNeighboors2(list, node, was_already):
    for x in net.adj[node]:
        if x.name:
            if not x.name in was_already:
                list.append(x.name.split('_2020')[0].split('2020_')[-1])
                was_already.append(x.name)
        else:
            for y in net.adj[x]:
                if y.name:
                    if not y.name in was_already:
                        list.append(y.name.split('_2020')[0].split('2020_')[-1])
                        was_already.append(y.name)
                else:
                    for z in net.adj[y]:
                        if z.name:
                            if not z.name in was_already:
                                list.append(z.name.split('_2020')[0].split('2020_')[-1])
                                was_already.append(z.name)



    return list, node,was_already
dict = {}
key_list = []
val_list = []

was_already = []
for x in net.nodes:

    if x.name:
        list = []
        was_already.append(x.name)
        FindNeighboors2(list,x,was_already)
        dict[x.name.split('_2020')[0].split('2020_')[-1]] = list
zeros = 0
testing=[]
for key in dict:
    for x in dict[key]:
        testing.append(x)


# from ete3 import Tree
#
# t = Tree('data/v.nh')
#
# relatives = {}
#
#
# def find_relatives(depth=2):
#     relatives = {}
#     for node in t.traverse("postorder"):
#         if node.name:
#             sub_rel = {}
#             u = node
#             for i in range(depth):
#                 u = u.up
#                 if u != None:
#                     sub_rel[i] = []
#                     for leaf in u:
#                         sub_rel[i].append(leaf.name)
#
#             relatives[node.name] = sub_rel
#     return relatives
#
#
# relatives = find_relatives(depth=2)

childless = []
for key in dict:
    if len(dict[key]) == 0:
        childless.append(key)

for key in childless:
    del dict[key]

names = subprocess.check_output("cat Cleaned_up_genes.fasta | grep EPI_ISL_", shell=True).split()

names = [name.decode("UTF-8") for name in names]

def get_number(header):
    return re.search("EPI_ISL_\\d+", header).group()

names = {get_number(name): name for name in names}

dict = {names[get_number(key)] : [names[get_number(el)] for el in dict[key]] for key in dict}

count = 0
for key in dict:
    count += len(dict[key])

file = open("extra_comparisons.txt", "w")
file.write(str(count) + "\n")
for parent in dict:
    for child in dict[parent]:
        file.write(parent)
        file.write("\n")
        file.write(child)
        file.write("\n")
file.close()