using LaplacianOpt
using JuMP
using CPLEX
using LinearAlgebra

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(

 
"solution_type" => "exact",
"optimizer" => "cplex",
"presolve" => true,
"optimizer_log" => true, 
"relax_integrality" => false,
                            
)

#------------------------------#
#      Optimization model      #
#------------------------------#
lom_optimizer = get_solver(params)
# data = LaplacianOpt.get_data(params)


# model_qc  = LaplacianOpt.build_LOModel(data)
# result_qc = LaplacianOpt.optimize_LOModel!(model_qc, optimizer = qcm_optimizer)
# LaplacianOpt.visualize_LOModel_solution(result_qc, data)
