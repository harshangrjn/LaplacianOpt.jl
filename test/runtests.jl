using LaplacianOpt
using Memento
using JuMP
using LinearAlgebra
using Test
using GLPK

const LO = LaplacianOpt
const LA = LinearAlgebra

# Suppress warnings during testing
LO.logger_config!("error")

@testset "LaplacianOpt" begin
    
    include("utility_tests.jl")
    include("lo_model_tests.jl")

end