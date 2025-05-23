import LaplacianOpt as LOpt
using JuMP
using Gurobi
# using GLPK

include("optimizers.jl")

#-------------------------------------#
#            MIP solver               #
# (> Gurobi 12.0.2 performs the best) #
#-------------------------------------#
lopt_optimizer = get_gurobi(solver_log = false)

#-------------------------------------#
#      User-defined input graphs      #
#-------------------------------------#
#= 
 Option I: Let the package parse a JSON file and obtain the data dictionary (data_dict)
 Sample JSON format: Check examples/instances
=#

# Spanning tree dataset
function data_I(num_nodes::Int, instance::Int)
    # Data format has to be as given in this JSON file
    file_path =
        joinpath(@__DIR__, "instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)
    augment_budget = (num_nodes - 1) # spanning tree constraint
    return data_dict, augment_budget
end

# Loading SLAM dataset
function data_SLAM(instance::String)
    # Data format has to be as given in this JSON file
    file_path = joinpath(@__DIR__, "instances/SLAM_dataset/$(instance).json")
    data_dict = LOpt.parse_file(file_path)
    augment_budget = 10
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
    augment_budget = 2
    return data_dict, augment_budget
end

#-------------------------------#
#      User-defined params      #
#-------------------------------#
num_nodes = 8
instance = 1
data_dict, augment_budget = data_I(num_nodes, instance)
# data_dict, augment_budget = data_SLAM("CSAIL") # "CSAIL", "intel", "ais2klinik"
# data_dict, augment_budget = data_II()

params = Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

#----------------------------------------------------------------#
#      Optimization model and visualize solution (optional)      #
#   Graph plots can be located inside `examples/plots` folder    #
#----------------------------------------------------------------#

# For more model options, check https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/src/types.jl
model_options = Dict{Symbol,Any}(
    :eigen_cuts_sizes => [data_dict["num_nodes"], 2],
    :topology_flow_cuts => true,
    :solution_type => "optimal",
)

result = LOpt.run_LOpt(
    params,
    lopt_optimizer;
    options = model_options,
    visualize_solution = false,  # Make it true to plot the graph solution
    visualizing_tool = "tikz",   # "tikz", "graphviz" are the options
)
