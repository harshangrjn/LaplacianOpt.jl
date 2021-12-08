import LaplacianOpt as LOpt
using JuMP
using CPLEX
using Gurobi
#using GLPK

include("solver.jl")

result_all = zeros(10,3)
lopt_optimizer = get_gurobi()

for i=1:1
    #-------------------------------#
    #      User-defined inputs      #
    #-------------------------------#
    params = Dict{String, Any}(
        "num_nodes" => 8,
        "instance" => i,
        "eigen_cuts_full" => true,
        "soc_linearized_cuts" => true
        )

    #-------------------------------------------#
    #      Optimization model and solution      #
    #-------------------------------------------#
    result = LOpt.run_LOpt_model(params, lopt_optimizer)
    
    result_all[i,1] = i 
    result_all[i,2] = result["objective"]
    result_all[i,3] = result["solve_time"]
end