# COACH by by Min Wu, Xiaoli Li, Chee-Keong Kwoh, and See-Kiong Ng. 
# "A core-attachment based method to detect protein complexes in ppi networks."

# Python 2 Author: True Price <jtprice@cs.unc.edu>
# Python 3 protclus version: Paul Scherer

from collections import defaultdict
from itertools import combinations
from functools import reduce
from py27hash.dict import Dict
from py27hash.set import Set
from tqdm import tqdm

from .cluster_alg import ClusterAlg

# return average degree and density for a graph
def graph_stats(graph):
    avg_deg = sum(len(n) for n in graph.values()) / float(len(graph))
    density = avg_deg / (len(graph)-1)
    return avg_deg, density

# return core nodes, given a graph and its average degree
get_core_nodes = lambda g,avg: Set(v for v,n in g.items() if len(n) >= avg)

# return NA score
NA_score = lambda a,b: float(len(a & b)**2) / (len(a) * len(b))


def core_removal(threshold, graph):
    if len(graph) == 1: # need at least two nodes in the graph...
        return [graph]

    avg_deg, density = graph_stats(graph)
    if density >= threshold:
        return [graph]
    else:
        # find and remove core nodes; create connected subcomponents
        core_nodes = get_core_nodes(graph, avg_deg)
        result = []
        subgraphs = []
        for v,n in graph.items():
            if v in core_nodes: continue
            n = n - core_nodes # note that we're reassigning n
            for s in subgraphs:
                if not n.isdisjoint(s):
                    s |= n
                    break
            else:
                subgraphs.append(n | Set([v]))
        # connected subcomponent joining
        i = 0
        while i < len(subgraphs) - 1:
            j = i + 1
            while j < len(subgraphs):
                if not subgraphs[i].isdisjoint(subgraphs[j]):
                    subgraphs[i] |= subgraphs[j]
                    subgraphs.pop(j)
                else:
                    j += 1
            i += 1
        # recursive core removal
        for s in subgraphs:
            tresults = core_removal(threshold, Dict((v,graph[v] & s) for v in s))
            for tc in tresults:
                nodes = Set()
                for v,n in tc.items():
                    nodes.add(v)
                    n |= graph[v] & core_nodes
                for c in core_nodes:
                    tc[c] = graph[c] & (nodes | core_nodes)
            result += tresults
        return result

class COACH(ClusterAlg):
    """Class for running and administrating the COACH algorithm"""
    def __init__(self, filename, density_threshold=0.7, affinity_threshold=0.225, closeness_threshold=0.5):
        super(COACH, self).__init__(filename)
        self.density_threshold = density_threshold
        self.affinity_threshold = affinity_threshold
        self.closeness_threshold = closeness_threshold

    def cluster(self, verbose=False):
        
        data = Dict()

        with open(self.filename, 'r') as f:
            for line in f:
                a,b = line.split()[:2]

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

        # step 1: find preliminary cores
        SC = [] # currently-detected preliminary cores
        count = 0
        for vertex,neighbors in tqdm(data.items()):
            # build neighborhood graph
            vertices = Set([vertex]) | neighbors
            size1_neighbors = Set()
            graph = { }
            for v in vertices:
                n = data[v] & vertices
                if len(n) > 1: # ignore size-1 vertices
                    graph[v] = n
                else:
                    size1_neighbors.add(v)
            if len(graph) < 2: # not enough connections in this graph
                continue
            graph[vertex] -= size1_neighbors

            # get core graph
            avg_deg,density = graph_stats(graph)
            core_nodes = get_core_nodes(graph, avg_deg)
            vertices = Set(graph.keys())
            for v in vertices - core_nodes:
                del graph[v]
            for n in graph.values():
                n &= core_nodes
            if len(graph) < 2: # not enough connections in this graph
                continue
            graph_nodes = Set(graph)

            # inner loop
            for sg in core_removal(self.density_threshold, graph):
                while True:
                    _,density = graph_stats(sg)
                    # if density threshold met, stop; else, remove min degree node
                    if density >= self.density_threshold: break
                    w = min(sg.items(), key=lambda k: len(k[1]))[0]
                    del sg[w]
                    for n in sg.values():
                        n.discard(w)

                sg_nodes = Set(sg)
                while graph_nodes - sg_nodes:
                    w = max(graph_nodes - sg_nodes,
                            key=lambda v: len(graph[v] & sg_nodes))
                    new_sg = sg.copy()
                    for v,n in new_sg.items():
                        if w in graph[v]:
                            n.add(w)
                    new_sg[w] = graph[w] & sg_nodes
                    _,density = graph_stats(new_sg)
                    if density < self.density_threshold: break
                    sg = new_sg
                    sg_nodes.add(w)

                # redundancy filtering
                max_sim = -1
                for i in range(len(SC)):
                    sim = NA_score(Set(SC[i]), sg_nodes)
                    if sim > max_sim:
                        max_sim = sim
                        index = i
                if max_sim < self.affinity_threshold:
                    SC.append(sg)
                else:
                    _,density_i = graph_stats(SC[index])
                    if density * len(sg) > density_i * len(SC[index]):
                        SC[index] = sg

        # step 2: adding peripheral proteins
        clusters = Set()
        for core in SC:
            nodes = frozenset(core)
            neighbors = reduce(lambda x,y: x|y, (data[v] for v in nodes)) - nodes
            neighbors -= Set(v for v in neighbors
              if float(len(data[v] & nodes)) / len(nodes) <= self.closeness_threshold)
            clusters.add(tuple(nodes | neighbors))


        self.clusters = clusters

        print ("Found %d clusters/protein complexes" % (len(clusters)))
        return clusters

# if __name__ == '__main__':
#     filename = "../data/unweighted_example_network.txt"
#     c = COACH(filename)
#     c.cluster()