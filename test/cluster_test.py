import numpy as np
import networkx as nx
from protclus import MCODE, DPCLUS

unweighted_filename = "data/unweighted_example_network.txt"


def test_mcode_cluster():
    """
    ## Testing MCODE on unweighted network
    """
    c = MCODE(unweighted_filename)
    assert len(c.clusters) == 0
    c.cluster()
    assert len(c.clusters) == 59
    assert len(c.clusters[-1]) == 4

def test_dpclus_cluster():
    """
    ## Testing DPCLUS on unweighted network
    """
    c = DPCLUS(unweighted_filename)
    assert len(c.clusters)==0
    c.cluster()
    assert len(c.clusters) == 901
    assert len(c.clusters[-1]) == 2