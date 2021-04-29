using LaplacianOpt
using JuMP
using CPLEX
using LinearAlgebra

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
"num_nodes" => 8,
"instance" => 1,
"data_type" => "old",

"eigen_cuts_full" => true,
 
"solution_type" => "exact",
"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true, 
"relax_integrality" => false
                        
)

#------------------------------#
#      Optimization model      #
#------------------------------#
lom_optimizer = get_solver(params)
data = LaplacianOpt.get_data(params)


model_lo  = LaplacianOpt.build_LOModel(data)
result_lo = LaplacianOpt.optimize_LOModel!(model_lo, optimizer = lom_optimizer)
# LaplacianOpt.visualize_LOModel_solution(result_qc, data)
