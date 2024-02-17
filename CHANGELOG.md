LaplacianOpt.jl Change Log
=========================

### Staged
- Solver logging option added in `optimizer.jl`
- Added option `minors_on_augment_edges` to include principal minor cuts only corresponding to vertices with augmentable edges - helps in reducing run times
- Cleaned up `get_minor_idx` function to make it more efficient 
- Function `_PMinorIdx` gets updated accordingly

### v0.6.2 
- Minor fix for testing eigen cut orthogonality if `projected_eigen_cuts` is active

### v0.6.1
- Readme update for dark/light theme logo display
- Fixed Docs compile issue
- Updated 100 and 60 nodes instances with solutions
- Minor updates in graphviz logging in `log.jl`
- `plot_solutions.jl` script added to plot graphs

### v0.6.0
- Added mutiple heuristics to handle both spanning trees and graphs with loops 
- Refactored `log.jl` to handle solutions from heuristics
- Included more user options for heuristic in `model_options`
- Cleaned up populating and logging of eigen cuts
- Updated docs and unit tests to reflect above changes

### v0.5.0
- Generalized eigen cuts to handle any size of non-negative minors (`2x2`-`NxN`)
- Dropped support for `constraint_eigen_cuts_on_2minors` and `constraint_eigen_cuts_on_3minors`
- Cleaned up project eigen-cuts, and applies to all sub-matrix sizes
- Added more dense-graph instances in `src/examples` folder (60 nodes)
- Added `src/examples/solutions` which conntains optimal/best-known solutions for all the dense-graph instances
- Updated docs and unit tests to reflect above changes

### v0.4.1
- Added more dense-graph instances, including the 9, 25 and 40 nodes instances in `src/examples` folder. 
- Added reference to the paper https://arxiv.org/abs/2304.08571

### v0.4.0
- Moved all the model options (including cuts) to `types.jl` into `LaplacianOptModelOptions`, a struct form 
- Streamlined default options for `LaplacianOptModelOptions`
- Clean up in `data.jl`
- Added evaluation of cheeger constant
- Added option to input `best_incumbent` to warm-start the solver
- Added certificate of optimality for MISDP problem
- Bug fix in cycle elimination cuts and clean up of flow cuts
- Added Pkg dependency
- Updated docs and unit tests to reflect above changes

### v0.3.2
- Logo update in README.md
- README.md added in the `src/examples`

### v0.3.1
- Minor change in `src/variables.jl`
- Added JuliaFormatter.toml and formatting workfow

### v0.3.0
- Includes adjacency of base and augments graphs in results dictionary
- Constraints added to support cycle graphs with max algberaic connectivity using `hamiltonian_cycle` in `graph_type`
- `_is_flow_cut_valid` can handle variable number of edges to be verified in cutset
- Added support for subtour elimination constraints
- Added option for `time_limit` in params (default = 10800 s)
- Clean up in `data.jl` for logging

### v0.2.1
- Update in `log.jl` for handling close to zero integral solutions before rounding
- Update in `lopt_model.jl` to handle displaying edge wts in plotting 
- Updated function references (for docs) in `data.jl` and `utility.jl`

### v0.2.0
- Breaking changes: Major restructuring of data format and generalization for max. algebraic connectivity edge augmentation problem 
- Data now takes a base graph with existing edges and the edges which can be augmented
- User input on the budget of edges to be augmented 
- data.jl updates to support new, generic data format 
- Bug fix and clean-up in topology flow cuts for connected components
- Added support for eigen cuts on 2x2 minors
- Added support for eigen cuts on 3x3 minors
- Added logo to the package
- Docs and tests updated to reflect above changes

### v0.1.8
- Added support for MOI v1.0+ 
- Added support for Graphs v1.0+
- Added support for JuMP v1.0+

### v0.1.7
- Added support for Graphs v1.5+, 1.6+
- Added support for JuMP v1.0+
- Added support for MOI v1.1+

### v0.1.6
- Fixed Tikzgraph issue in tests 
- Transitioned from LightGraphs to Graphs
- Updated and cleaned up `examples/optimizer.jl`
- `src/data.jl` cleanup for handling `"optimizer"`
- Minor docs update 

### v0.1.5
- Added CITATION.bib
- Added support for Gurobi MIP solver in `examples/optimizer.jl` 
- Updated `LO` to `LOpt`
- `lopt_model.jl` function calls updated with `LOpt`
- Minor docs cleanups

### v0.1.4 
- Updated types to handle graph structs, `GraphData`
- Updated docs for utility functions
- Updated units tests to cover `GraphData`

### v0.1.3
- Types bug fix for constructor in `LaplacianOptModel`
- Updated tests in max span tree model 
- Updated cleaned-up docs

### v0.1.2
- Clean up of plot functions in log.jl
- Clean up of docs 
- Bug fix in tests for plots

### v0.1.1
- Added utility functions for evaluation of algebraic connectivity 
- Added utility functions for evaluation of fieder vector
- Support for SOC relaxations for upper bounds 
- Documentation stable update
- More unit tests for SOC and utility functions

### v0.1.0
- Initial function-based working implementation 
- Added preliminary documentation 
- Added preliminary unit tests
- Github Actions set-up
