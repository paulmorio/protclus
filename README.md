# protclus - Minimal Python 3 library for Protein Complex Discovery in PPI Networks

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![PyPI](https://img.shields.io/pypi/v/protclus)
![PyPI - License](https://img.shields.io/pypi/l/protclus)
[![Build Status](https://travis-ci.org/paulmorio/protclus.svg?branch=master)](https://travis-ci.org/paulmorio/protclus)
![Codecov](https://img.shields.io/codecov/c/github/paulmorio/protclus)

This is a Python 3 library containing clustering algorithms chiefly used for protein complex discovery in protein-protein interaction (PPI) networks. 

It is inspired by the collection of Python 2 scripts in https://github.com/trueprice/python-graph-clustering (which is unfortunately due to deprecate in 2021) --- but updated for Python 3 and rewritten to as a set of classes and methods. As a result it is significantly easier to import and integrate into projects. Additionally it benefits from some additional flexibility in the output and some improvements on the speed of the algorithms.

It benefits from extremely simple usage that is shared between all of the methods

```python
# Example using MCODE
from protclus import MCODE
filename = "data/unweighted_example_network.txt"
c = MCODE(filename)
c.cluster()

print(c.clusters)
c.save_clusters("mcode_example_clusters.txt")
```

## Methods

### Currently Included

- ***MCODE*** by Gary D. Bader and Christopher W. V. Hogue "An automated method for finding molecular complexes in large protein interaction networks."
- ***DPCLUS*** by Md Altaf-Ul-Amin et al. "Development and implementation of an algorithm for detection of protein complexes in large interaction networks."
- ***IPCA*** by Min Li, Jian-er Chen, Jian-xin Wang, Bin Hu, and Gang Chen. Modifying the dpclus algorithm for identifying protein complexes based on new topological structures.
- ***COACH*** by by Min Wu, Xiaoli Li, Chee-Keong Kwoh, and See-Kiong Ng. "A core-attachment based method to detect protein complexes in ppi networks."

### Coming Soon
- ***Graph Entropy*** by E. C. Kenley and Y. Cho "Entropy-Based Graph Clustering: Application to Biological and Social Networks"

## Installation
Installation from PyPI

```bash
pip install protclus
```

Installation from source from the project root folder where `setup.py` can be done via
```bash
pip install -e .
```

### Testing

```bash
python setup.py test
```

## Citation
If this work was of use to you please cite the original authors of each of the algorithms and the previous Python2 script authors.

Additionally please consider citing the following paper, as this library was developed as part of

```
@misc{protclus,
      title={Incorporating network based protein complex discovery into automated model construction}, 
      author={Paul Scherer and Maja Trȩbacz and Nikola Simidjievski and Zohreh Shams and Helena Andres Terre and Pietro Liò and Mateja Jamnik},
      year={2020},
      eprint={2010.00387},
      archivePrefix={arXiv},
      primaryClass={q-bio.MN}
}
```
