@testset "Algebraic Connectivity: Optimal solution tests" begin

    num_nodes = 5
    instance  = 1
    file_path = joinpath(@__DIR__, "..","examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)

    params = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => (num_nodes-1),
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

    for i=1:num_nodes
        L_val[i,i] += result_lo["solution"]["γ_var"] * (num_nodes - 1)/num_nodes
    end

    for i=1:(num_nodes - 1)
        for j=(i+1):num_nodes
            L_val[i,j] -= result_lo["solution"]["γ_var"] / num_nodes
            L_val[j,i] = L_val[i,j]
        end
    end

    @test isapprox(LA.eigvals(L_val)[2], result_lo["objective"])
    @test isapprox(LA.eigvals(L_val)[2], result_lo["solution"]["γ_var"])
end

@testset "Max Span Tree: Optimal solution tests" begin
    
    num_nodes = 5
    instance  = 1
    file_path = joinpath(@__DIR__, "..","examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)

    params = Dict{String, Any}(
        "data_dict"      => data_dict,
        "augment_budget" => (num_nodes-1)
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

    num_nodes = 15
    instance  = 1
    file_path = joinpath(@__DIR__, "..","examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)

    params = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => (num_nodes-1)
        )

    result_mst = LaplacianOpt.run_MaxSpanTree(params, glpk_optimizer, lazy_callback=true)
    
    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 3625.89701005, atol=1E-6)
    @test isapprox(result_mst["solution"]["z_var"][1,7], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][2,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][3,4], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][4,14], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][4,15], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][5,14], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][6,15], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][7,12], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][7,15], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][8,13], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][8,14], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][9,10], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][9,13], 1.0)
    @test isapprox(result_mst["solution"]["z_var"][11,12], 1.0)

    # Test multi commodity flow
    result_mst_1 = LaplacianOpt.run_MaxSpanTree(params, glpk_optimizer, lazy_callback=false)
    @test result_mst_1["termination_status"] == MOI.OPTIMAL
    @test result_mst_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst_1["objective"], 3625.89701005, atol=1E-6)

end

@testset "SOC relaxations - 1: constraint_soc_cuts_on_2minors" begin

    num_nodes = 5
    instance  = 5
    file_path = joinpath(@__DIR__, "..","examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)

    params_1 = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => (num_nodes-1),
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
        "data_dict" => data_dict,
        "augment_budget" => (num_nodes-1),
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

@testset "Minor-based eigen relaxations: 2x2 and 3x3" begin

    num_nodes = 5
    instance  = 3
    file_path = joinpath(@__DIR__, "..","examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
    data_dict = LOpt.parse_file(file_path)

    params_1 = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => (num_nodes-1),
        "eigen_cuts_full" => true,
        "eigen_cuts_2minors"  => true,
        "eigen_cuts_3minors"  => true
    )

    result_1 = LaplacianOpt.run_LOpt(params_1, glpk_optimizer)
    
    @test result_1["termination_status"] == MOI.OPTIMAL
    @test result_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_1["objective"], 8.356739455, atol=1E-6)
    @test isapprox(result_1["solution"]["z_var"][1,4], 1.0)
    @test isapprox(result_1["solution"]["z_var"][2,4], 1.0)
    @test isapprox(result_1["solution"]["z_var"][2,5], 1.0)
    @test isapprox(result_1["solution"]["z_var"][3,4], 1.0)
end

@testset "Test base graph with existing edges" begin
    function data_I()
        data_dict = Dict{String, Any}()
        data_dict["num_nodes"] = 4
        data_dict["adjacency_base_graph"] = [0 2 0 0
                                             2 0 3 0 
                                             0 3 0 4
                                             0 0 4 0]
        data_dict["adjacency_augment_graph"] = [0 0 4 8
                                                0 0 0 7 
                                                4 0 0 0
                                                8 7 0 0]
        augment_budget = 2
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_I()

    params = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => augment_budget,
        "eigen_cuts_full" => true,
        "soc_linearized_cuts" => true,
        "eigen_cuts_2minors"  => true,
        "eigen_cuts_3minors"  => true
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer)
    
    @test result["termination_status"] == MOI.OPTIMAL
    @test result["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result["objective"], 7.8551986404, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1,4], 1.0)
    @test isapprox(result["solution"]["z_var"][2,4], 1.0)

    function data_II()
        num_nodes = 5
        instance  = 1
        file_path = joinpath(@__DIR__, "..", "examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance)_test.json")
        data_dict = LOpt.parse_file(file_path)
        augment_budget = 3
        return data_dict, augment_budget
    end

    data_dict, augment_budget = data_II()

    params = Dict{String, Any}(
        "data_dict" => data_dict,
        "augment_budget" => augment_budget
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer)
    
    @test result["termination_status"] == MOI.OPTIMAL
    @test result["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result["objective"], 34.404789725, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1,4], 1.0)
    @test isapprox(result["solution"]["z_var"][2,4], 1.0)
    @test isapprox(result["solution"]["z_var"][2,5], 1.0)
end