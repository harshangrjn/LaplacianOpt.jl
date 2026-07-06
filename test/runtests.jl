using LaplacianOpt
using Test

import JuMP
import LinearAlgebra
import Logging
import GLPK
import MathOptInterface

const LOpt = LaplacianOpt
const LA = LinearAlgebra
const MOI = MathOptInterface

# Suppress warnings during testing
LOpt.logger_config!("error")

function _test_warn(f, msg)
    io = IOBuffer()
    old_logger = LOpt._LOGGER[]
    LOpt._LOGGER[] = Logging.ConsoleLogger(io, Logging.Warn)
    try
        f()
    finally
        LOpt._LOGGER[] = old_logger
    end

    return @test occursin(msg, String(take!(io)))
end

function _test_throws_with_logger(f, error_type)
    io = IOBuffer()
    old_logger = LOpt._LOGGER[]
    LOpt._LOGGER[] = Logging.ConsoleLogger(io, Logging.Error)
    try
        @test_throws error_type f()
    finally
        LOpt._LOGGER[] = old_logger
    end
end

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

test_time_limit() = 100

@testset "LaplacianOpt" begin
    @testset "logging helpers" begin
        _test_warn(
            () -> LOpt._PMinorIdx(3, [4], false, Dict{String,Any}()),
            "Detected maximum eigen-cut size",
        )
        _test_throws_with_logger(
            () -> LOpt._catch_data_input_error(2, 3, 1, 1.0),
            ErrorException,
        )
        LOpt.logger_config!("error")
    end

    include("utility_tests.jl")
    include("lo_model_tests.jl")
    include("heuristics_tests.jl")
    include("relaxations_test.jl")
    include("log_tests.jl")
end
