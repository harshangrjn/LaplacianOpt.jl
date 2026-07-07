<h1 align="center" margin=0px>
  <!-- <img src="https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/docs/src/assets/logo_header_light.svg#gh-light-mode-only" width=75%>
  <img src="https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/docs/src/assets/logo_header_dark.svg#gh-dark-mode-only"   width=75%>
  <br> -->
  <a href="https://github.com#gh-light-mode-only">
  <img src="https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/docs/src/assets/logo_header_light.svg#gh-light-mode-only" width=75%>
  </a>
  <a href="https://github.com#gh-dark-mode-only">
    <img src="https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/docs/src/assets/logo_header_dark.svg#gh-dark-mode-only" width=75%>
  </a>
  <br>
  A Julia Package for Maximizing Algebraic Connectivity of Graphs
</h1>

Status: 
[![CI](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl/branch/main/graph/badge.svg?token=7EKATOHLYL)](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl)
[![Documentation](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/LaplacianOpt.jl/dev/)
[![version](https://juliahub.com/docs/LaplacianOpt/version.svg)](https://juliahub.com/ui/Packages/LaplacianOpt/V1JEg/)
## LaplacianOpt
**LaplacianOpt** is a Julia package for designing robust weighted graphs by maximizing the algebraic connectivity of their Laplacian matrices. Given a weighted base graph (possibly empty), a set of candidate weighted edges, and an augmentation budget `K`, LaplacianOpt selects `K` edges that make the resulting graph as well connected as possible.

[Algebraic connectivity](https://dml.cz/bitstream/handle/10338.dmlcz/101168/CzechMathJ_23-1973-2_11.pdf) is the second smallest eigenvalue of the graph Laplacian. Larger values indicate stronger global connectivity, and this measure is widely used to study robustness, synchronizability, and sparsification of networks.

The package supports mixed-integer formulations, polyhedral relaxations, cutting-plane methods, and heuristic approaches for finding globally optimal and high-quality locally optimal solutions. A common use case is designing sparse but robust measurement networks for robot localization under a limited edge budget; starting from `N` vertices, no existing edges, a complete candidate graph, and `K = N - 1`, LaplacianOpt searches for a spanning tree with maximum algebraic connectivity.

## Usage
Clone the repository and start Julia in the project environment:

```bash
julia --project=.
```

From the Julia package prompt, instantiate dependencies and run the test suite:

```julia
pkg> instantiate
pkg> test
```

See the `examples` directory for optimization scripts, input data, plotting utilities, and sample workflows.

## Bug reports and support
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/LaplacianOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc.

## Acknowledgement
This work was supported by Los Alamos National Laboratory (LANL)'s LDRD Early Career Research Award (20190590ECR) and [LANL-TAMU's collaborative research project](https://nationallabsoffice.tamus.edu/the-texas-am-university-system-and-los-alamos-national-laboratory-partner-to-design-robust-networks/) grant. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citing LaplacianOpt
If you find LaplacianOpt.jl useful in your work, we request you to cite the following papers [\[link-1\]](https://doi.org/10.1109/ECC.2015.7330770) [\[link-2\]](https://doi.org/10.1109/TCNS.2024.3431408): 
```bibtex
@article{LOpt_TCNS2024,
  title={Optimal robust network design: Formulations and algorithms for maximizing algebraic connectivity},
  author={Somisetty, Neelkamal and Nagarajan, Harsha and Darbha, Swaroop},
  journal={IEEE Transactions on Control of Network Systems},
  url = {https://doi.org/10.1109/TCNS.2024.3431408},
  pages={918--929},
  volume={12},
  number={1},
  year={2024},
  publisher={IEEE}
}

@inproceedings{LOpt_ECC2015,
  title={On maximizing algebraic connectivity of networks for various engineering applications},
  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},
  booktitle={European Control Conference (ECC)},
  pages={1626--1632},
  year={2015},
  organization={IEEE}
}
```