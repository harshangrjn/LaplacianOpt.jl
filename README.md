<p align="center">
<img width="790px" src="https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/logo.png" alt="https://github.com/harshangrjn/LaplacianOpt.jl/tree/master/docs/src/assets/docs_header_dark.png"/>
</p>

Status: 
[![CI](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl/branch/main/graph/badge.svg?token=7EKATOHLYL)](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl)
[![Documentation](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/LaplacianOpt.jl/dev/)
## LaplacianOpt
**LaplacianOpt** is a Julia package which implements polyhedral relaxation-based algorithms for the maximimum algebraic connectivity augmentation problem on weighted graph Laplacians. More specifically, given a weighted base graph with existing edges (could be empty), a set of candidate weighted edges for augmentation, and an augmentation budget (`K`), this package finds a set of `K` edges to augment to the base graph such that the resulting graph has maximum algebraic conenctivity with optimality guarantees. For example, if the base graph with `N` vertices is empty, set of candidate edges form a complete graph, and `K = (N-1)`, this packages finds a spanning tree with maximum algebraic connectivity.

[Algebraic connectivity](https://dml.cz/bitstream/handle/10338.dmlcz/101168/CzechMathJ_23-1973-2_11.pdf) is the second smallest eigenvalue of the graph Laplacian. The magnitude of this value reflects how well connected the overall graph is. This connectivity measure has been used in analyzing the robustness and synchronizability of complex networks, and in graph sparsification techniques. 

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

Check the "examples" folder on how to use this package.

## Bug reports and support
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/LaplacianOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc.

## Acknowledgement
This work was supported by Los Alamos National Laboratory (LANL)'s LDRD Early Career Research Award (20190590ECR) and [LANL-TAMU's collaborative research project](https://nationallabsoffice.tamus.edu/the-texas-am-university-system-and-los-alamos-national-laboratory-partner-to-design-robust-networks/) grant. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citing LaplacianOpt
If you find LaplacianOpt.jl useful in your work, we request you to cite the following papers [\[link-1\]](https://doi.org/10.1109/ECC.2015.7330770) [\[link-2\]](https://doi.org/10.1115/1.4028955): 
```bibtex
@inproceedings{LOpt_ECC2015,
  title={On maximizing algebraic connectivity of networks for various engineering applications},
  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},
  booktitle={European Control Conference (ECC)},
  pages={1626--1632},
  year={2015},
  organization={IEEE}
}

@article{LOpt_ASME2015,
  title={Synthesizing robust communication networks for unmanned aerial vehicles with resource constraints},
  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},
  journal={Journal of Dynamic Systems, Measurement, and Control},
  volume={137},
  number={6},
  pages={061001},
  year={2015},
  publisher={American Society of Mechanical Engineers}
}
```