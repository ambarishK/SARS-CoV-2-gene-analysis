from Bio import Phylo

import networkx, pylab

tree = Phylo.read('test.nh', 'newick')

net = Phylo.to_networkx(tree)

list = networkx.generate_adjlist(net)
#print(net.clade_relations)
print(list)
counter = 0

import sys
sys.setrecursionlimit(150000)
def FindNeighboors(list, node):
    for x in net.adj[node]:
        if x.name:
            list.append(x.name)
        else:
            FindNeighboors(list,x)
    return list, node

dict = {}
key_list = []
val_list = []


for x in net.nodes:
    list = []
    FindNeighboors(list,x)
    break
    if  x.name:

        list = []

        #retarded version


        
        # for y in net.adj[x]:
        #     for z in net.adj[y]:
        #         if z.name  and z.name != x.name:
        #             counter += 1
        #             list.append((z.name.split('_2020')[0].split('2020_')[-1]))
        #             val_list.append(z.name)
        #
        #         else:
        #             for w in net.adj[z]:
        #                 if w.name and w.name !=x.name:
        #                     list.append((w.name.split('_2020')[0].split('2020_')[-1]))
        #                     val_list.append(w.name)




        dict[x.name.split('_2020')[0].split('2020_')[-1]] = list
        key_list.append(x.name)
zeros = 0

print(dict)
