# IPCA by Min Li, Jian-er Chen, Jian-xin Wang, Bin Hu, and Gang Chen. 
# Modifying the dpclus algorithm for identifying protein complexes 
# based on new topological structures.

# Python 2 Author: True Price <jtprice@cs.unc.edu>
# Python 3 protclus version: Paul Scherer

import sys
from itertools import combinations
from collections import defaultdict
from py27hash.dict import Dict
from py27hash.set import Set

from .cluster_alg import ClusterAlg

class IPCA(ClusterAlg):
    """
    Class for running and administrating the IPCA clustering algorithm

    """

    def __init__(self, filename, t_in=0.5):
        super(IPCA, self).__init__(filename)
        self.t_in = t_in

    def cluster(self, verbose=False):
        # data = defaultdict(Set) # node id => neighboring node ids

        data = Dict()
        # read in graph
        with open(self.filename, 'r') as f:
            counter = 0
            for line in f:
                a,b = line.split()[:2]
                counter += 1
                if a in data:
                    data[a].add(b)
                else:
                    data[a] = Set()
                    data[a].add(b)
                if b in data:
                    data[b].add(a)
                else:
                    data[b] = Set()
                    data[b].add(a)

        # weights = defaultdict(int)
        weights = Dict()
        for a,b in combinations(data, 2):
            if b not in data[a]: continue   
            shared = len(data[a] & data[b])
            if a in weights:
                weights[a] += shared
            else:
                weights[a] = 0
                weights[a] += shared
            if b in weights:
                weights[b] += shared
            else:
                weights[b] = 0
                weights[b] += shared

        unvisited = Set(data)
        num_clusters = 0
        clusters = []
        
        # print(unvisited)
        # return 0

        # Potential culprit
        seed_nodes = sorted(data, key=lambda k: (weights[k],len(data[k])), reverse=True)

        for seed in seed_nodes: # get highest degree node
            if seed not in unvisited: continue

            cluster = Set((seed,next(iter(data[seed])))) # seed and random neighbor

            while True:
                # rank neighbors by the number of edges between the node and cluster nodes
                frontier = sorted((len(data[p] & cluster),p) for p in
                    Set.union(*((data[n] - cluster) for n in cluster)))

                # do this until IN_vk < T_IN, SP <= 2 is met, or no frontier nodes left
                found = False
                while frontier and not found:
                    m_vk,p = frontier.pop()
                    if m_vk < self.t_in * len(cluster): break
                    c_2neighbors = data[p] & cluster
                    c_2neighbors.update(*(data[c] & cluster for c in c_2neighbors))
                    if cluster == c_2neighbors:
                        found = True
                        break

                if not found: break
                    
                # otherwise, add the node to the cluster
                cluster.add(p)

            unvisited -= cluster

            if verbose:
                print (' '.join(cluster))

            num_clusters += 1

            if verbose:
                print (num_clusters, len(cluster), len(unvisited))

            clusters.append(cluster)

            if not unvisited: break

        self.clusters = clusters

# if __name__ == '__main__':
#     filename = "../data/unweighted_example_network.txt"
#     c = IPCA(filename)
#     c.cluster()