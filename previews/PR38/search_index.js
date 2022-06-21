var documenterSearchIndex = {"docs":
[{"location":"function_references/#LaplacianOpt.jl-Function-References","page":"Function References","title":"LaplacianOpt.jl Function References","text":"","category":"section"},{"location":"function_references/","page":"Function References","title":"Function References","text":"Modules = [LaplacianOpt]","category":"page"},{"location":"function_references/#LaplacianOpt.GraphData","page":"Function References","title":"LaplacianOpt.GraphData","text":"GraphData\n\nThe composite mutable struct, GraphData, holds matrices of adjacency matrix, laplacian matrix, fiedler vector and algebraic connectivity.\n\n\n\n\n\n","category":"type"},{"location":"function_references/#LaplacianOpt.LaplacianOptModel","page":"Function References","title":"LaplacianOpt.LaplacianOptModel","text":"LaplacianOptModel\n\nThe composite mutable struct, LaplacianOptModel, holds dictionaries for input data, abstract JuMP model for optimization, variable references and result from solving the JuMP model.\n\n\n\n\n\n","category":"type"},{"location":"function_references/#LaplacianOpt._violated_eigen_vector-Tuple{Array{Float64}}","page":"Function References","title":"LaplacianOpt._violated_eigen_vector","text":"_violated_eigen_vector(W::Array{Float64}; tol = 1E-6)\n\nGiven a square symmetric matrix, this function returns an eigen vector corresponding to it's most  violated eigen value w.r.t positive semi-definiteness of the input matrix.  \n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.algebraic_connectivity-Tuple{Array{Float64}}","page":"Function References","title":"LaplacianOpt.algebraic_connectivity","text":"algebraic_connectivity(adjacency_matrix::Array{Float64})\n\nReturns the algebraic connectivity or the   second smallest eigenvalue of the Laplacian matrix, for an input weighted adjacency matrix of the graph. \n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.fiedler_vector-Tuple{Array{Float64}}","page":"Function References","title":"LaplacianOpt.fiedler_vector","text":"fiedler_vector(adjacency_matrix::Array{Float64})\n\nReturns the Fiedler vector or the eigenvector corresponding to the  second smallest eigenvalue of the Laplacian matrix for an input weighted adjacency matrix of the graph. \n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.get_objective_bound-Tuple{JuMP.Model}","page":"Function References","title":"LaplacianOpt.get_objective_bound","text":"\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.get_objective_value-Tuple{JuMP.Model}","page":"Function References","title":"LaplacianOpt.get_objective_value","text":"\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.laplacian_matrix-Tuple{Array{Float64}}","page":"Function References","title":"LaplacianOpt.laplacian_matrix","text":"laplacian_matrix(adjacency_matrix::Array{Float64})\n\nReturns the weighted Laplacian matrix for an input weighted adjacency matrix of the graph. \n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.logger_config!-Tuple{Any}","page":"Function References","title":"LaplacianOpt.logger_config!","text":"allows the user to set the logging level without the need to add Memento\n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.optimal_graph_edges-Tuple{Array{Float64}}","page":"Function References","title":"LaplacianOpt.optimal_graph_edges","text":"optimal_graph_edges(adjacency_matrix::Array{Float64})\n\nReturns a vector of tuples of edges corresponding  to an input adjacency matrix of the graph. \n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.relaxation_bilinear-Tuple{JuMP.Model, JuMP.VariableRef, JuMP.VariableRef, JuMP.VariableRef}","page":"Function References","title":"LaplacianOpt.relaxation_bilinear","text":"relaxation_bilinear(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef)\n\nGeneral relaxation of a binlinear term using McCormick relaxations. This can be used to obtain specific variants  in partiuclar cases of variables (when either or both the variables are binary)\n\nz >= JuMP.lower_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.lower_bound(x)*JuMP.lower_bound(y)\nz >= JuMP.upper_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.upper_bound(x)*JuMP.upper_bound(y)\nz <= JuMP.lower_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.lower_bound(x)*JuMP.upper_bound(y)\nz <= JuMP.upper_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.upper_bound(x)*JuMP.lower_bound(y)\n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.silence-Tuple{}","page":"Function References","title":"LaplacianOpt.silence","text":"Suppresses information and warning messages output by LaplacianOpt, for fine grained control use the Memento package\n\n\n\n\n\n","category":"method"},{"location":"function_references/#LaplacianOpt.variable_domain-Tuple{JuMP.VariableRef}","page":"Function References","title":"LaplacianOpt.variable_domain","text":"variable_domain(var::JuMP.VariableRef)\n\nComputes the valid domain of a given JuMP variable taking into account bounds and the varaible's implicit bounds (e.g. binary).\n\n\n\n\n\n","category":"method"},{"location":"quickguide/#Quick-Start-Guide","page":"Quick Start guide","title":"Quick Start Guide","text":"","category":"section"},{"location":"quickguide/#Getting-started","page":"Quick Start guide","title":"Getting started","text":"","category":"section"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"After the installation of LaplacianOpt and CPLEX (use GLPK if open-source is preferable) from the Julia package manager, and the cost matrix of the complete graph (for example, see examples/instances/5_nodes) is provided in the JSON format, an optimization model to maximize the algebraic connectivity of the weighted graph's Laplacian matrix can be executed with a few lines of code as follows:","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"import LaplacianOpt as LOpt\nusing JuMP\nusing Gurobi\n\nparams = Dict{String, Any}(\n\"num_nodes\" => 5,\n\"instance\" => 1\n)\n\nlopt_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, \"presolve\" => 1) \nresults = LOpt.run_LOpt(params, lopt_optimizer)","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"tip: Tip\nRun times of LaplacianOpt's mathematical optimization models are significantly faster using Gurobi as the underlying mixed-integer programming (MIP) solver. Note that this solver's individual-usage license is available free for academic purposes. ","category":"page"},{"location":"quickguide/#Extracting-results","page":"Quick Start guide","title":"Extracting results","text":"","category":"section"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"The run commands (for example, run_LOpt) in LaplacianOpt return detailed results in the form of a dictionary. This dictionary can be used for further processing of the results. For example, for the given instance of a complete graph, the algorithm's runtime and the optimal objective value (maximum algebraic connectivity) can be accessed with,","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"results[\"solve_time\"]\nresults[\"objective\"]","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"The \"solution\" field contains detailed information about the solution produced by the optimization model. For example, one can obtain the edges of the optimal graph toplogy from the symmetric adjacency matrix with,","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"optimal_graph = LOpt.optimal_graph_edges(results[\"solution\"][\"z_var\"])","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"Further, algebraic connectivity and the Fiedler vector of the optimal graph topology (though can be applied on any graph topology) can be obtained from the adjacency matrix with,","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"optimal_adjacency = results[\"solution\"][\"z_var\"] .* LOpt.get_data(params)[\"edge_weights\"] \ngraph_data = LOpt.GraphData(optimal_adjacency)\nprintln(\"Algebraic connectivity: \", graph_data.ac)\nprintln(\"Fiedler vector: \", graph_data.fiedler)","category":"page"},{"location":"quickguide/#Visualizing-results","page":"Quick Start guide","title":"Visualizing results","text":"","category":"section"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"LaplacianOpt also currently supports the visualization of optimal graphs layouts obtained from the results dictionary (from above. To do so, these are the two options: ","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"TikzGraphs package for a simple and quick visualization of the graph layout without support to include edge weights, which can be executed with ","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"data = LOpt.get_data(params)\nLOpt.visualize_solution(results, data, visualizing_tool = \"tikz\")","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"Graphviz package for better visualization of weighted graphs. To this end, LaplacianOpt generates the raw .dot file, which can be further visualized using the Graphviz software either via the direct installation on the computer or using an online front-end visualization GUI (for example, see Edotor). Dot files can be generated in LaplacianOpt with ","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"data = LOpt.get_data(params)\nLOpt.visualize_solution(results, data, visualizing_tool = \"graphviz\")","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"For example, on a weighted complete graph with 10 nodes in instance #1, the optimal spanning tree with maximum algebraic connectivity, out of 10^8 feasible solutions, obtained by LaplacianOpt (using Graphviz visualization) is shown below ","category":"page"},{"location":"quickguide/","page":"Quick Start guide","title":"Quick Start guide","text":"(Image: Optimal solution)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"<align=\"center\"/>\n<img width=\"790px\" class=\"display-light-only\" src=\"assets/docs_header.png\" alt=\"assets/docs_header.png\"/>\n<img width=\"790px\" class=\"display-dark-only\" src=\"assets/docs_header_dark.png\" alt=\"assets/docs_header.png\"/>","category":"page"},{"location":"#Documentation","page":"Introduction","title":"Documentation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CurrentModule = LaplacianOpt","category":"page"},{"location":"#Overview","page":"Introduction","title":"Overview","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"LaplacianOpt is a Julia package which implements polyhedral relaxation-based algorithms for optimization of weighted graph Laplacians. Given a complete weighted, undirected graph, this package provides an optimal (and if preferable, an approximate) spanning tree which has the maximum second smallest eigenvalue of the graph Laplacian, also known as the Algebraic Connectivity of the graph. This package also implements various types of convex relaxations to the Laplacian optimization problem. ","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To use LaplacianOpt, first download and install Julia. Note that the current version of LaplacianOpt is compatible with Julia 1.0 and later. ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"The latest stable release of LaplacianOpt can be installed using the Julia package manager with","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(\"LaplacianOpt\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"At least one mixed-integer programming solver is required for running LaplacianOpt. The well-known CPLEX or the Gurobi solver is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale graphs. However, the open-source GLPK solver is also compatible with LaplacianOpt which can be installed via the package manager with","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(\"GLPK\")","category":"page"},{"location":"#Unit-Tests","page":"Introduction","title":"Unit Tests","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"To run the tests in the package, run the following command after installing the LaplacianOpt package.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.test(\"LaplacianOpt\")","category":"page"},{"location":"#Citing-LaplacianOpt","page":"Introduction","title":"Citing LaplacianOpt","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"If you find LaplacianOpt.jl useful in your work, we request you to cite the following papers [link-1] [link-2]: ","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"@inproceedings{LOpt_ECC2015,\n  title={On maximizing algebraic connectivity of networks for various engineering applications},\n  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},\n  booktitle={European Control Conference (ECC)},\n  pages={1626--1632},\n  year={2015},\n  organization={IEEE}\n}\n\n@article{LOpt_ASME2015,\n  title={Synthesizing robust communication networks for unmanned aerial vehicles with resource constraints},\n  author={Nagarajan, Harsha and Rathinam, Sivakumar and Darbha, Swaroop},\n  journal={Journal of Dynamic Systems, Measurement, and Control},\n  volume={137},\n  number={6},\n  pages={061001},\n  year={2015},\n  publisher={American Society of Mechanical Engineers}\n}","category":"page"}]
}
