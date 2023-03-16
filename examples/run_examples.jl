import LaplacianOpt as LOpt
using JuMP
using CPLEX
# using Gurobi
# using GLPK

include("optimizer.jl")

#----------------------------------#
#            MIP solver            #
# (Cplex 22.1 performs the best)   #
#----------------------------------#
lopt_optimizer = get_cplex()

#-------------------------------------#
#      User-defined input graphs      #
#-------------------------------------#
#= 
 Option I: Let the package parse a JSON file and obtain the data dictionary (data_dict)
 Sample JSON format: Check examples/instances
=#
function data_I(num_nodes::Int, instance::Int)
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
    augment_budget = 2
    return data_dict, augment_budget
end

#-------------------------------#
#      User-defined params      #
#-------------------------------#
gamma_lb_8 = [22.804157, 24.320666, 26.411124, 28.691222, 22.505054, 25.216707, 22.875174, 28.439690, 26.796540, 27.491292]
gamma_lb_10 = [34.237096, 41.448828, 37.730852, 41.461819, 34.319331, 39.972654, 36.165099, 42.329072, 39.403384, 34.916098]
gamma_lb_12 = [54.0522, 53.2107, 47.2228, 43.9330, 51.1286, 56.9622, 57.2901, 53.2338, 53.5628, 50.6987]

result_tab = zeros(10, 3)
num_nodes = 8

for k = 1:10
    println("Instance number: ", k)
    instance = k
    data_dict, augment_budget = data_I(num_nodes, instance)

    params = Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    #----------------------------------------------------------------#
    #      Optimization model and visualize solution (optional)      #
    #   Graph plots can be located inside `examples/plots` folder    #
    #----------------------------------------------------------------#
   
    # For more model options, check https://github.com/harshangrjn/LaplacianOpt.jl/blob/master/src/types.jl
    model_options = Dict{Symbol,Any}(
        :eigen_cuts_full => false,
        :eigen_cuts_2minors => true,
        :eigen_cuts_2minors => true,
        :topology_flow_cuts => true,
        :cheeger_cuts => false,
        :best_lower_bound => eval(Symbol(:gamma_lb_, num_nodes))[k],
        :solution_type => "optimal",
    )

    result = LOpt.run_LOpt(
        params,
        lopt_optimizer;
        options = model_options,
        visualize_solution = false,  # Make this true to plot the graph solution
        visualizing_tool = "tikz",   # "graphviz" is another option
    )

    result_tab[k,1] = k
    result_tab[k,2] = result["objective"]
    result_tab[k,3] = (result["solve_time"] - result["cheeger_cuts_separation_time"])
    @show result_tab[k,:]
end
