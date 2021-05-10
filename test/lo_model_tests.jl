@testset "Algebraic Connectivity: Optimal solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/solver.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "data_type" => "old",
         
        "solution_type" => "exact",
        "optimizer" => "glpk",
        "presolve" => true,
        "optimizer_log" => true, 
        "relax_integrality" => false                      
    )

    lom_optimizer = get_solver(params)
    result_lo = LaplacianOpt.run_LOpt_model(params, lom_optimizer)

    @test result_lo["termination_status"] == MOI.OPTIMAL
    @test result_lo["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_lo["objective"], 11.1592, atol=1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1,4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4,1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2,4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4,2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3,4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4,3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3,5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5,3], 1.0)

    @test LA.issymmetric(result_lo["solution"]["W_var"])
    
    L_val = result_lo["solution"]["W_var"]

    for i=1:params["num_nodes"]
        L_val[i,i] += result_lo["solution"]["γ_var"] * (params["num_nodes"] - 1)/params["num_nodes"]
    end

    for i=1:(params["num_nodes"] - 1)
        for j=(i+1):params["num_nodes"]
            L_val[i,j] -= result_lo["solution"]["γ_var"] / params["num_nodes"]
            L_val[j,i] = L_val[i,j]
        end
    end

    @test isapprox(LA.eigvals(L_val)[2], result_lo["objective"])
    @test isapprox(LA.eigvals(L_val)[2], result_lo["solution"]["γ_var"])
end

@testset "Max Span Tree: Optimal solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/solver.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "optimizer" => "glpk"
        )

    lom_optimizer = get_solver(params)
    result_mst = LaplacianOpt.run_MaxSpanTree_model(params, lom_optimizer)

    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 113.7588, atol=1E-4)
    @test isapprox(result_mst["solution"]["z_var"][1,2], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][2,1], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][2,5], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][5,2], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][3,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][4,3], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][3,5], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][5,3], 1.0)

end
