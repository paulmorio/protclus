# DPCLUS algorithm from :Altaf-Ul-Amin M, Shinbo Y, Mihara K, Kurokawa K, Kanaya S: 
# Development and implementation of an algorithm for detection of protein complexes 
# in large interaction networks. BMC Bioinformatics 2006

# Python 2 Author: True Price <jtprice@cs.unc.edu>
# Python 3 protclus version: Paul Scherer

import sys
from itertools import combinations
from collections import defaultdict

D_THRESHOLD = 0.9
CP_THRESHOLD = 0.5

# dictionary type that returns zero for missing values
# used here in 'edges' dictionary
class zerodict(dict):
    def __missing__(self, k):
        return 0

def dpclus(filename):
    data = defaultdict(set) # node id => neighboring node ids

    # to match the original DPClus output, we need to keep track of the order in
    # which nodes were read from the file, in order to maintain cluster order when
    # breaking ties in our priority ranking (values are decreasing, i.e. negative)
    node_index = dict()
    b_list = [] # the original algorithm indexes column b after column a... >:-(

    # read in graph
    with open(filename, 'r') as f:
        for line in f:
            a,b = line.split()[:2]
            data[a].add(b)
            data[b].add(a)
            node_index.setdefault(a, -len(node_index))
            #node_index.setdefault(b, -len(node_index))
            b_list.append(b)
    for b in b_list:
        node_index.setdefault(b, -len(node_index))

    unvisited = set(data)
    num_clusters = 0

    clusters = []
    while unvisited:
        # get highest degree node
        seed = max(unvisited, key=lambda k: (len(data[k]&unvisited),node_index[k]))
        frontier = data[seed] & unvisited
        if not frontier: break # no connections left to analyze

        edges,weights = defaultdict(zerodict), defaultdict(int)
        for a,b in combinations(unvisited, 2):
            if b not in data[a]: continue 
            shared = len(data[a] & data[b] & unvisited)
            edges[a][b],edges[b][a] = shared, shared
            weights[a] += shared
            weights[b] += shared

        max_w,_,node = max((w,node_index[n],n) for n,w in weights.items())
        if max_w > 0:
            seed = node
            frontier = data[seed] & unvisited

        cluster = set((seed,))
        cluster_degrees = {seed: 0}
        nn,ne = 1, 0 # number of nodes, edges in cluster

        print(seed, end=" ")

        while frontier:
            # find higest priority node:
            # 0. number of edges between node and cluster nodes
            # 1. sum of edge weights between node and cluster nodes
            # 2. the node's index
            # 3. the node itself
            e_nk,_,_,p = max((len(data[n] & cluster), sum(edges[n][c] for c in cluster), node_index[n], n) for n in frontier)

            density = 2. * (ne + e_nk) / (nn * (nn+1))
            if density < D_THRESHOLD:
                break # adding the node gives too low density; cluster is finished

            # if all nodes only have one connecting edge in cluster, use "fine-tuning"
            # this orders by the number of neighbors in the frontier, minus the
            # connectedness of the attached (cluster) node within the cluster
            cp = CP_THRESHOLD
            if e_nk == 1 and len(cluster) > 1:
                print("::")
                n_degree = dict() # node => fine-tuning parameter
                for n in frontier:
                    for c in cluster: # find adjacent cluster node
                        if n in edges[c]: break
                    n_degree[n] = len(data[n] & frontier) - cluster_degrees[c]
                p = max(frontier, key=lambda k: (n_degree[k],node_index[k]))
                if n_degree[p] > 0:
                    cp /= 2.
            if (e_nk / density / (nn+1)) < cp:
                break # no good node found; cluster is finished

            print(p, end=" ")

            # otherwise, add the node to the cluster
            cluster.add(p)
            nn,ne = (nn+1), (ne+e_nk)

            cluster_degrees[p] = e_nk
            for n in data[p] & cluster:
                cluster_degrees[n] += 1

            frontier = set.union(*((data[n] - cluster) & unvisited for n in cluster))

        # add overlapping nodes
        # frontier[2] stores our fine-tuning parameter, in this case
        frontier_nodes = (set.union(*(data[c] for c in cluster)) - cluster)
        frontier = sorted([len(data[n] & cluster), sum(edges[n][c] for c in cluster), 0, node_index[n], n] for n in frontier_nodes)
        
        # "fine-tuning"
        fine_tuning = False
        if frontier and frontier[-1][0] == 1:
            for n in frontier:
                for c in cluster: # find adjacent cluster node
                    if n[4] in edges[c]: break
                n[2] = len(data[n[4]] & frontier_nodes) - cluster_degrees[c]
            frontier.sort(key=lambda k: (k[2],k[3]))
            fine_tuning = True
        
        # iterate through visited neighbors and update accordingly
        while frontier:
            e_nk,_,w,_,p = frontier.pop()
            cp = CP_THRESHOLD  / 2. if fine_tuning and w > 0 else CP_THRESHOLD

            density = 2. * (ne + e_nk) / (nn * (nn+1))
            if density < D_THRESHOLD or (e_nk / density / (nn+1)) < cp: continue

            print (p, end=" ")

            # add node to the cluster
            cluster.add(p)
            nn,ne = (nn+1), (ne+e_nk)

            cluster_degrees[p] = e_nk
            for n in data[p] & cluster:
                cluster_degrees[n] += 1

            # update E_nk for the other nodes on the frontier
            for n in frontier:
                if p in edges[n[4]]:
                    n[0] += 1

        unvisited -= cluster

        num_clusters += 1
        print (num_clusters, nn, 2. * ne / nn / (nn-1))

        clusters.append(cluster)

if __name__ == '__main__':
    import sys
    dpclus(sys.argv[1])
 
