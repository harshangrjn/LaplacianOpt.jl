import LaplacianOpt as LOpt
using JuMP
using CPLEX
using Gurobi
#using GLPK

include("solver.jl")

result_all = zeros(10,3)
for i=1:10
    #-------------------------------#
    #      User-defined inputs      #
    #-------------------------------#
    params = Dict{String, Any}(
        "num_nodes" => 8,
        "instance" => i,
        "eigen_cuts_full" => true,
        "soc_linearized_cuts" => true,
        "optimizer" => "gurobi"
        )

    #-------------------------------------------#
    #      Optimization model and solution      #
    #-------------------------------------------#
    lom_optimizer = get_solver(params)
    result = LOpt.run_LOpt_model(params, lom_optimizer)
    
    result_all[i,1] = i 
    result_all[i,2] = result["objective"]
    result_all[i,3] = result["solve_time"]
end