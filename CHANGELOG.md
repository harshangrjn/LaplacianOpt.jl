LaplacianOpt.jl Change Log
=========================


### v0.1.6
- Fixed Tikzgraph issue for tests 
- Transitioned from LightGraphs to Graphs
- Updated and cleaned up `examples/solver.jl`
- `src/data.jl` cleanup for handling `"optimizer"`
- Minor docs update

### v0.1.5
- Added CITATION.bib
- Added support for Gurobi MIP solver in `examples/solver.jl` 
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
