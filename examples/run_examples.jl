import LaplacianOpt as LOpt
using JuMP
using CPLEX
# using Gurobi
# using GLPK

include("optimizer.jl")

#---------------------------------------------#
#         User-defined MIP solver             #
# (Cplex 20.1/22.1 seems to perform the best) #
#---------------------------------------------#
lopt_optimizer = get_cplex()

#-------------------------------------#
#      User-defined input graphs      #
#-------------------------------------#
#= 
 Option I: Let the package parse a JSON file and obtain the data dictionary (data_dict)
 Sample JSON format: Check examples/instances
=#
function data_I()
    num_nodes = 8
    instance = 1
    # Data format has to be as given in this JSON file
    file_path =
        joinpath(@__DIR__, "instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)
    augment_budget = (num_nodes - 1) # spanning tree constraint
    return data_dict, augment_budget
end

#=  
 Option II: Directly input the data dictionary (data_dict) with 
            num_nodes, adjacency_base_graph and adjacency_augment_graph
=#
function data_II()
    data_dict = Dict{String,Any}()
    data_dict["num_nodes"] = 4
    data_dict["adjacency_base_graph"] = [0 2 0 0; 2 0 3 0; 0 3 0 4; 0 0 4 0]
    data_dict["adjacency_augment_graph"] = [0 0 4 8; 0 0 0 7; 4 0 0 0; 8 7 0 0]
    augment_budget = 3
    return data_dict, augment_budget
end

#-------------------------------#
#      User-defined params      #
#-------------------------------#
data_dict, augment_budget = data_I()

params = Dict{String,Any}(
    "data_dict" => data_dict,
    "augment_budget" => augment_budget,
    "eigen_cuts_full" => true,
    "soc_linearized_cuts" => false,
    "eigen_cuts_2minors" => false,
    "eigen_cuts_3minors" => false,
    "topology_flow_cuts" => true,
    # "time_limit"          => 3600,
)

#----------------------------------------------------------------#
#      Optimization model and visualize solution (optional)      #
#   Graph plots can be located inside `examples/plots` folder    #
#----------------------------------------------------------------#
result = LOpt.run_LOpt(
    params,
    lopt_optimizer;
    visualize_solution = false,  # Make this true to plot the graph solution
    visualizing_tool = "tikz", # "graphviz" is another option
    display_edge_weights = false,
)
