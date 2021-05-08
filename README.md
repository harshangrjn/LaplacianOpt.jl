Status: 
[![CI](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml)
[![Documentation](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/documentation.yml)
[![codecov](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl/branch/main/graph/badge.svg?token=7EKATOHLYL)](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl)
## LaplacianOpt
**LaplacianOpt.jl** is a Julia package which implements polyhedral relaxation-based algorithms for optimization of weighted graph Laplacians. Given a complete weighted graph, this package provides an optimal (and an approximate) spanning tree which has the maximum second largest eigenvalue of the graph Laplacian, also known as the algebraic connectivity of the graph. This package also implements various types of relaxations to the Laplacian optimization problem. 

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

Check the "examples" folder on how to use this package.

## Bug reports and support
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/LaplacianOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc.

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research Award, *"20190590ECR: Discrete Optimization Algorithms for Provable Optimal Quantum Circuit Design"*. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citation
If you find LaplacianOpt.jl useful in your work, we request you to cite the following paper [\[link\]](https://doi.org/10.1109/ECC.2015.7330770): 
```bibtex
@inproceedings{NagarajanRathinamDarbha2015,
  title={On maximizing algebraic connectivity of networks for various engineering applications},
  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},
  booktitle={European Control Conference (ECC)},
  pages={1626--1632},
  year={2015},
  organization={IEEE}
}
```