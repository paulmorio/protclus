# General clustering algorithm class
# Defines some behaviour common to all of the clustering algorithms

# Author: Paul Scherer
# MIT LICENSE

class ClusterAlg(object):
    """General class for clustering algorithms defines some features common to all 
    clustering algorithms in the package

    """

    def __init__(self, filename):
        self.filename = filename
        self.clusters = []

    def __str__(self):
        return (f"Clustering algorithm on {self.filename}")

    def save_clusters(self, filehandle):
        """Saves clusters, one cluster per line into the input filehandle"""
        with open(filehandle, 'w') as fh:
            for c in self.clusters:
                fh.write(' '.join(c) + "\n")
