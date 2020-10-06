# Bader et al. (2003) MCODE clustering algorithm from
# "An automated method for finding molecular complexes in large protein interaction networks"

# Python 2 Author: True Price <jtprice@cs.unc.edu>
# Python 3 protclus version: Paul Scherer

import sys
from collections import defaultdict
from tqdm import tqdm

from cluster_alg import ClusterAlg


class MCODE(ClusterAlg):
    """Class for running and administrating Bader et al.'s MCODE algorithm"""

    def __init__(self, filename, weight_threshold=0.2):
        super(MCODE, self).__init__(filename)
        self.weight_threshold = 1 - weight_threshold

    def cluster(self):
        edges = defaultdict(set)  # node id => neighboring node ids

        # Read edgelist
        with open(self.filename, 'r') as f:
            for line in f:
                a, b = line.split()[:2]
                edges[a].add(b)
                edges[b].add(a)
        print ('## Input graph loaded; %i nodes' % (len(edges),))

        # Clusters list
        clusters = []

        # Stage 1: Vertex Weighting
        print ('## Weighting vertices...')
        weights = dict((v, 1.) for v in edges)
        for i, v in tqdm(enumerate(edges)):
            neighborhood = set((v,)) | edges[v]
            # if node has only one neighbor, we know everything we need to know
            if len(neighborhood) <= 2:
                continue

            # see if larger k-cores exist
            k = 1  # highest valid k-core
            while neighborhood:
                k_core = neighborhood.copy()
                invalid_nodes = True
                while invalid_nodes and neighborhood:
                    invalid_nodes = set(n for n in neighborhood if len(
                        edges[n] & neighborhood) <= k)
                    neighborhood -= invalid_nodes
                k += 1  # on exit, k will be one greater than we want

            # vertex weight = k-core number * density of k-core
            weights[v] = (k - 1) * (sum(len(edges[n] & k_core)
                                        for n in k_core) / (2. * len(k_core)**2))

        # Stage 2: Molecular Complex Prediction
        print('## Molecular complex prediction...')
        unvisited = set(edges)
        num_clusters = 0

        for seed in sorted(weights, key=weights.get, reverse=True):
            if seed not in unvisited:
                continue

            cluster, frontier = set((seed,)), set((seed,))
            w = weights[seed] * self.weight_threshold
            while frontier:
                cluster.update(frontier)
                unvisited -= frontier
                frontier = set(n for n in set.union(
                    *(edges[n] for n in frontier)) & unvisited if weights[n] > w)

            # Haircut: only keep 2-core complexes
            invalid_nodes = True
            while invalid_nodes and cluster:
                invalid_nodes = set(
                    n for n in cluster if len(edges[n] & cluster) < 2)
                cluster -= invalid_nodes

            if cluster:
                # fluff never really seems to improve anything...
                # cluster.update(
                # n for n in set.union(*(edges[c] for c in cluster)) & unvisited
                # if densities[n] > FLUFF_THRESHOLD)

                print (' '.join(cluster))
                num_clusters += 1
                print (num_clusters, len(cluster), seed)
                clusters.append(cluster)

        self.clusters = clusters

# if __name__ == '__main__':
#     filename = "../data/unweighted_example_network.txt"
#     c = MCODE(filename)
#     c.cluster()
#     c.save_clusters("/tmp/mcode_test.txt")
