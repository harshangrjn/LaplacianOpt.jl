@testset "Algebraic Connectivity: Optimal solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/optimizer.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "data_type" => "old",
         
        "solution_type" => "exact",
        "presolve" => true,
        "optimizer_log" => true, 
        "relax_integrality" => false                      
    )

    result_lo = LaplacianOpt.run_LOpt(params, glpk_optimizer)

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
    
    file_path = joinpath(@__DIR__,"..", "examples/optimizer.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1
        )

    result_mst = LaplacianOpt.run_MaxSpanTree(params, glpk_optimizer)

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

@testset "Max Span Tree: Lazy callback tests" begin

    file_path = joinpath(@__DIR__,"..", "examples/optimizer.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 15,
        "instance" => 5,
        "data_type" => "old",
        "eigen_cuts_full" => false,
        "topology_flow_cuts" => true
        )

    result_mst = LaplacianOpt.run_MaxSpanTree(params, glpk_optimizer, lazy_callback=true)
    
    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 3576.56202331, atol=1E-6)
    @test isapprox(result_mst["solution"]["z_var"][1,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][2,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][3,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][4,15], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][7,15], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][8,15], 1.0)
end

@testset "SOC relaxations - 1: constraint_soc_cuts_on_2minors" begin

    file_path = joinpath(@__DIR__,"..", "examples/optimizer.jl")
    include(file_path)

    params_1 = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 5,
        "eigen_cuts_full" => false,
        "soc_linearized_cuts" => true,
        "topology_flow_cuts" => true
    )

    result_1 = LaplacianOpt.run_LOpt(params_1, glpk_optimizer)
    
    @test result_1["termination_status"] == MOI.OPTIMAL
    @test result_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_1["objective"], 15.992792245, atol=1E-6)
    @test isapprox(result_1["solution"]["z_var"][1,5], 1.0)
    @test isapprox(result_1["solution"]["z_var"][2,5], 1.0)
    @test isapprox(result_1["solution"]["z_var"][3,4], 1.0)
    @test isapprox(result_1["solution"]["z_var"][3,5], 1.0)
    @test isapprox(result_1["solution"]["z_var"][4,3], 1.0)

    params_2 = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 5,
        "eigen_cuts_full" => true,
        "soc_linearized_cuts" => true
    )

    result_2 = LaplacianOpt.run_LOpt(params_2, glpk_optimizer)
    
    @test result_2["termination_status"] == MOI.OPTIMAL
    @test result_2["primal_status"] == MOI.FEASIBLE_POINT
    @test result_2["objective"] <= result_1["objective"]
    @test isapprox(result_2["objective"], 13.065372790, atol=1E-6)
    @test isapprox(result_1["solution"]["z_var"], result_2["solution"]["z_var"], atol = 1E-6)
end