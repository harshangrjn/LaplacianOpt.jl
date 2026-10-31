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
[![codecov](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl/branch/master/graph/badge.svg?token=7EKATOHLYL)](https://codecov.io/gh/harshangrjn/LaplacianOpt.jl)
[![Documentation](https://github.com/harshangrjn/LaplacianOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/LaplacianOpt.jl/dev/)
[![version](https://img.shields.io/github/v/release/harshangrjn/LaplacianOpt.jl?label=version)](https://juliahub.com/ui/Packages/LaplacianOpt/V1JEg/)
## LaplacianOpt
**LaplacianOpt** is a Julia package for designing robust weighted graphs by maximizing the algebraic connectivity of their Laplacian matrices. Given a weighted base graph (possibly empty), a set of candidate weighted edges, and an augmentation budget `K`, LaplacianOpt selects `K` edges that make the resulting graph as well connected as possible.

[Algebraic connectivity](https://dml.cz/bitstream/handle/10338.dmlcz/101168/CzechMathJ_23-1973-2_11.pdf) is the second smallest eigenvalue of the graph Laplacian. Larger values indicate stronger global connectivity, and this measure is widely used to study robustness, synchronizability, and sparsification of networks.

The package supports mixed-integer formulations, polyhedral relaxations, cutting-plane methods, and heuristic approaches for finding globally optimal and high-quality locally optimal solutions. A common use case is designing sparse but robust measurement networks for robot localization under a limited edge budget; starting from `N` vertices, no existing edges, a complete candidate graph, and `K = N - 1`, LaplacianOpt searches for a spanning tree with maximum algebraic connectivity.

## Graph sparsification for pose-graph SLAM

In pose-graph Simultaneous Localization and Mapping (SLAM) application, a limited loop-closure budget turns robot measurement selection into a graph-design problem. LaplacianOpt selects the loop closures that maximize the algebraic connectivity of the measurement graph, reinforcing weakly constrained regions and improving the accuracy of the final pose estimate (see Figure below). Because a larger algebraic connectivity lowers the Cramér–Rao bound on estimation error ([Boumal et al., 2014](https://doi.org/10.1093/imaiai/iat006)), a sparse graph with high algebraic connectivity can support accurate estimation for surface reconstruction. See the
[research documentation](https://harshangrjn.github.io/LaplacianOpt.jl/dev/research/#graph-sparsification-for-pose-graph-slam)
for the problem setup details.

<p align="center">
  <a href="https://harshangrjn.github.io/LaplacianOpt.jl/dev/research/#graph-sparsification-for-pose-graph-slam">
    <img src="docs/src/assets/posegraph_algebraic_connectivity.png" alt="Five pose-graph reconstructions showing improved recovery as algebraic connectivity increases" width="100%">
  </a>
</p>

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
Please report any issues via the GitHub **[issue tracker](https://github.com/harshangrjn/LaplacianOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc.

## Acknowledgement
This work was supported by Los Alamos National Laboratory (LANL)'s LDRD Early Career Research Award (20190590ECR) and [LANL-TAMU's collaborative research project](https://nationallabsoffice.tamus.edu/the-texas-am-university-system-and-los-alamos-national-laboratory-partner-to-design-robust-networks/) grant. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citing LaplacianOpt
If you find LaplacianOpt useful in your work, please cite the publication(s)
most relevant to your use of the package. The
[Research & Publications](https://harshangrjn.github.io/LaplacianOpt.jl/dev/research/)
page provides a chronological overview of the research, paper links, and
copyable BibTeX citations.