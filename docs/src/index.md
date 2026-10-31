# LaplacianOpt Documentation

```@meta
CurrentModule = LaplacianOpt
```
## Overview
**[LaplacianOpt](https://github.com/harshangrjn/LaplacianOpt.jl)** is a Julia package that implements polyhedral relaxation-based algorithms for the maximum algebraic connectivity augmentation problem on weighted graph Laplacians. More specifically, given a weighted base graph with existing edges (which could be empty), a set of candidate weighted edges for augmentation, and an augmentation budget (`K`), this package finds a set of `K` edges to augment the base graph such that the resulting graph has maximum algebraic connectivity with optimality guarantees. For example, given a base graph with `N` vertices and no edges, a set of candidate edges that forms a complete graph, and `K = N - 1`, this package finds a spanning tree with maximum algebraic connectivity.

[Algebraic connectivity](https://dml.cz/bitstream/handle/10338.dmlcz/101168/CzechMathJ_23-1973-2_11.pdf) is the second smallest eigenvalue of the graph Laplacian. The magnitude of this value reflects how well connected the overall graph is. This connectivity measure has been used in analyzing the robustness and synchronizability of complex networks, and in graph sparsification techniques. 

## Installation 
To use [LaplacianOpt](https://github.com/harshangrjn/LaplacianOpt.jl), first [download and install](https://julialang.org/downloads/) Julia. The current version of LaplacianOpt requires Julia 1.10 or later.

The latest stable release of [LaplacianOpt](https://github.com/harshangrjn/LaplacianOpt.jl) can be installed using the Julia package manager with

```julia
import Pkg
Pkg.add("LaplacianOpt")
```

At least one mixed-integer programming solver is required to run LaplacianOpt. [CPLEX](https://github.com/jump-dev/CPLEX.jl) and [Gurobi](https://github.com/jump-dev/Gurobi.jl) are recommended for solving larger graph instances. The open-source [GLPK](https://github.com/jump-dev/GLPK.jl) solver is also compatible with LaplacianOpt and can be installed through the package manager with

```julia
import Pkg
Pkg.add("GLPK")
```

## Unit Tests
To run the tests in the package, run the following command after installing the [LaplacianOpt](https://github.com/harshangrjn/LaplacianOpt.jl) package.

```julia
import Pkg
Pkg.test("LaplacianOpt")
```

## Citing LaplacianOpt
If you find [LaplacianOpt](https://github.com/harshangrjn/LaplacianOpt.jl) useful in your work, please cite the publication(s)
most relevant to your use of the package. The
[Research & Publications](research.md#Publications) page provides a chronological overview
of the research, paper links, and copyable BibTeX citations.