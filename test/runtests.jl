using LaplacianOpt
using Test

import Memento
import JuMP
import LinearAlgebra
import GLPK
import MathOptInterface

const LOpt = LaplacianOpt
const LA = LinearAlgebra
const MOI = MathOptInterface

# Suppress warnings during testing
LOpt.logger_config!("error")

glpk_optimizer = JuMP.optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true)

function data_spanning_tree(num_nodes::Int, instance::Int)
    num_nodes = num_nodes
    instance = instance
    file_path = joinpath(
        @__DIR__,
        "..",
        "examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json",
    )
    data_dict = LOpt.parse_file(file_path)
    augment_budget = num_nodes - 1
    return data_dict, augment_budget
end

test_time_limit() = 60

@testset "LaplacianOpt" begin
    include("utility_tests.jl")
    include("lo_model_tests.jl")
    include("heuristics_tests.jl")
    include("relaxations_test.jl")
    include("log_tests.jl")
end
