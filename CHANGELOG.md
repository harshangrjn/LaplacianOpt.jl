LaplacianOpt.jl Change Log
=========================

### v0.2.0
- Breaking changes: Major restructuring of data format and generalization for max. algebraic connectivity edge augmentation problem 
- Data now takes a base graph with existing edges and the edges which can be augmented
- User input on the budget of edges to be augmented 
- data.jl updates to support new, generic data format 
- Bug fix and clean-up in topology flow cuts for connected components
- Added logo to the package

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
