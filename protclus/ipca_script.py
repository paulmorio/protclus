# IPCA by Min Li, Jian-er Chen, Jian-xin Wang, Bin Hu, and Gang Chen. 
# Modifying the dpclus algorithm for identifying protein complexes 
# based on new topological structures.

# Python 2 Author: True Price <jtprice@cs.unc.edu>
# Python 3 protclus version: Paul Scherer


import sys
from itertools import combinations
from collections import defaultdict

T_IN = 0.5

# dictionary type that returns zero for missing values
# used here in 'edges' dictionary
class zerodict(dict):
    def __missing__(self, k):
        return 0

def ipca(filename):
    data = defaultdict(set) # node id => neighboring node ids

    # read in graph
    with open(filename, 'r') as f:
        for line in f:
            a,b = line.split()[:2]
            data[a].add(b)
            data[b].add(a)

    weights = defaultdict(int)
    for a,b in combinations(data, 2):
        if b not in data[a]: continue   
        shared = len(data[a] & data[b])
        weights[a] += shared
        weights[b] += shared

    unvisited = set(data)
    num_clusters = 0

    seed_nodes = sorted(data,key=lambda k: (weights[k],len(data[k])),reverse=True)

    for seed in seed_nodes: # get highest degree node
        if seed not in unvisited: continue

        cluster = set((seed,iter(data[seed]).next())) # seed and random neighbor

        while True:
            # rank neighbors by the number of edges between the node and cluster nodes
            frontier = sorted((len(data[p] & cluster),p) for p in
                set.union(*((data[n] - cluster) for n in cluster)))

            # do this until IN_vk < T_IN, SP <= 2 is met, or no frontier nodes left
            found = False
            while frontier and not found:
                m_vk,p = frontier.pop()
                if m_vk < T_IN * len(cluster): break
                c_2neighbors = data[p] & cluster
                c_2neighbors.update(*(data[c] & cluster for c in c_2neighbors))
                if cluster == c_2neighbors:
                    found = True
                    break

            if not found: break
                
            # otherwise, add the node to the cluster
            cluster.add(p)

        unvisited -= cluster
        print (' '.join(cluster))

        num_clusters += 1
        print (num_clusters, len(cluster), len(unvisited))

        if not unvisited: break

if __name__ == '__main__':
    import sys

    ipca(sys.argv[1])
 