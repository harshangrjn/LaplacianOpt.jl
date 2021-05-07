using LaplacianOpt
using Test

import Memento
import JuMP
import LinearAlgebra
import GLPK
import MathOptInterface

const LO = LaplacianOpt
const LA = LinearAlgebra
const MOI = MathOptInterface

# Suppress warnings during testing
LO.logger_config!("error")

@testset "LaplacianOpt" begin
    
    include("utility_tests.jl")
    include("lo_model_tests.jl")

end