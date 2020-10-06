# protclus - A small Python3 library for protein complex discovery with topological clustering techniques on PPI networks

This is a Python3 library containing clustering algorithms chiefly used for protein complex discovery in PPI networks. It is inspired by the collection of Python2 modules in https://github.com/trueprice/python-graph-clustering (which is unfortunately due to deprecate in 2021), but updated for Python3 and rewritten as a class-based library. As a result it is significantly easier to import and integrate into projects. Additionally it benefits from some additional flexibility in the output and some improvements on the speed of the algorithms. 

## Prerequisites
- Python 3.6+
- NetworkX 2.3+

## Installation
```bash
pip install protclus
```
## Citation
If this work was of use to you please cite the original authors of each of the algorithms and the previous Python2 module authors.
Additionally please consider citing the following paper, as this library was built in conjunction.

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
