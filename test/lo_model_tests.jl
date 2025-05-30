@testset "Algebraic Connectivity: Optimal solution tests" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [num_nodes],
        :projected_eigen_cuts => false,
        :topology_flow_cuts => true,
        :solution_type => "optimal",
        :relax_integrality => false,
        :time_limit => test_time_limit(),
    )
    result_lo = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result_lo["termination_status"] == MOI.OPTIMAL
    @test result_lo["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_lo["objective"], 11.1592, atol = 1E-4)

    for (i, j) in [(1, 4), (4, 1), (2, 4), (4, 2), (3, 4), (4, 3), (3, 5), (5, 3)]
        @test isapprox(result_lo["solution"]["z_var"][i, j], 1.0)
    end

    @test LA.issymmetric(result_lo["solution"]["W_var"])

    L_val = result_lo["solution"]["W_var"]

    for i in 1:num_nodes
        L_val[i, i] += result_lo["solution"]["γ_var"] * (num_nodes - 1) / num_nodes
    end

    for i in 1:(num_nodes-1)
        for j in (i+1):num_nodes
            L_val[i, j] -= result_lo["solution"]["γ_var"] / num_nodes
            L_val[j, i] = L_val[i, j]
        end
    end

    @test isapprox(LA.eigvals(L_val)[2], result_lo["objective"])
    @test isapprox(LA.eigvals(L_val)[2], result_lo["solution"]["γ_var"])
end

@testset "Max Span Tree: Optimal solution tests" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [],
        :formulation_type => "max_span_tree",
        :topology_flow_cuts => false,
        :topology_multi_commodity => true,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )
    result_mst = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 113.7588, atol = 1E-4)
    for (i, j) in [(1, 2), (2, 5), (3, 4), (3, 5)]
        @test isapprox(result_mst["solution"]["z_var"][i, j], 1.0)
        @test isapprox(result_mst["solution"]["z_var"][j, i], 1.0)
    end
end

@testset "Max Span Tree: Lazy callback tests" begin
    num_nodes = 15
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [],
        :formulation_type => "max_span_tree",
        :topology_flow_cuts => true,
        :topology_multi_commodity => false,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )
    result_mst = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 3625.89701005, atol = 1E-6)

    for (i, j) in [
        (1, 7),
        (2, 4),
        (3, 4),
        (4, 14),
        (4, 15),
        (5, 14),
        (6, 15),
        (7, 12),
        (7, 15),
        (8, 13),
        (8, 14),
        (9, 10),
        (9, 13),
        (11, 12),
    ]
        @test isapprox(result_mst["solution"]["z_var"][i, j], 1.0)
    end

    # Test multi commodity flow
    model_options[:topology_flow_cuts] = false
    model_options[:topology_multi_commodity] = true
    result_mst_1 = LaplacianOpt.run_LOpt(params, glpk_optimizer, options = model_options)
    @test result_mst_1["termination_status"] == MOI.OPTIMAL
    @test result_mst_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst_1["objective"], 3625.89701005, atol = 1E-6)
end

@testset "Max Span Tree: Flow cuts tests" begin
    num_nodes = 15
    instance = 3
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [],
        :formulation_type => "max_span_tree",
        :topology_flow_cuts => true,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )
    result_mst = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result_mst["termination_status"] == MOI.OPTIMAL
    @test result_mst["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 3695.02394015, atol = 1E-6)

    # Without flow cuts
    model_options[:max_cycle_length] = 0
    result_mst_1 = LaplacianOpt.run_LOpt(params, glpk_optimizer, options = model_options)
    @test result_mst_1["termination_status"] == MOI.OPTIMAL
    @test result_mst_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst["objective"], 3695.02394015, atol = 1E-6)

    # Without flow cuts
    model_options[:topology_flow_cuts] = false
    result_mst_2 = LaplacianOpt.run_LOpt(params, glpk_optimizer, options = model_options)
    @test result_mst_2["termination_status"] == MOI.OPTIMAL
    @test result_mst_2["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_mst_2["objective"], 3695.82875361, atol = 1E-6)
end

@testset "SOC relaxations - 1: constraint_soc_cuts_on_2minors" begin
    num_nodes = 5
    instance = 7
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params_1 =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [],
        :projected_eigen_cuts => false,
        :soc_linearized_cuts => true,
        :topology_flow_cuts => true,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result_1 = LaplacianOpt.run_LOpt(params_1, glpk_optimizer; options = model_options)

    @test result_1["termination_status"] == MOI.OPTIMAL
    @test result_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_1["objective"], 19.529478817, atol = 1E-6)
    for (i, j) in [(1, 2), (2, 5), (3, 4), (3, 5)]
        @test isapprox(result_1["solution"]["z_var"][i, j], 1.0, atol = 1E-6)
    end

    params_2 =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [num_nodes],
        :projected_eigen_cuts => false,
        :soc_linearized_cuts => true,
        :topology_flow_cuts => true,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result_2 = LaplacianOpt.run_LOpt(params_2, glpk_optimizer; options = model_options)

    @test result_2["termination_status"] == MOI.OPTIMAL
    @test result_2["primal_status"] == MOI.FEASIBLE_POINT
    @test result_2["objective"] <= result_1["objective"]
    @test isapprox(result_2["objective"], 12.41137427, atol = 1E-6)
end

@testset "Minor-based eigen relaxations: 2x2 and 3x3" begin
    num_nodes = 5
    instance = 3
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params_1 =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [2, 3, num_nodes],
        :projected_eigen_cuts => false,
        :topology_multi_commodity => true,
        :topology_flow_cuts => false,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result_1 = LaplacianOpt.run_LOpt(params_1, glpk_optimizer; options = model_options)

    @test result_1["termination_status"] == MOI.OPTIMAL
    @test result_1["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result_1["objective"], 8.356739455, atol = 1E-6)
    for (i, j) in [(1, 4), (2, 4), (2, 5), (3, 4)]
        @test isapprox(result_1["solution"]["z_var"][i, j], 1.0)
    end
end

@testset "Test base graph with existing edges" begin
    function data_II()
        data_dict = Dict{String,Any}()
        data_dict["num_nodes"] = 4
        data_dict["adjacency_base_graph"] = [
            0 2 0 0
            2 0 3 0
            0 3 0 4
            0 0 4 0
        ]
        data_dict["adjacency_augment_graph"] = [
            0 0 4 8
            0 0 0 7
            4 0 0 0
            8 7 0 0
        ]
        augment_budget = 2
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_II()

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [2, data_dict["num_nodes"]],
        :soc_linearized_cuts => true,
        :projected_eigen_cuts => false,
        :topology_flow_cuts => false,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["termination_status"] == MOI.OPTIMAL
    @test result["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result["objective"], 7.8551986404, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result["solution"]["z_var"][2, 4], 1.0)

    function data_II()
        num_nodes = 5
        instance = 1
        file_path = joinpath(
            @__DIR__,
            "..",
            "examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance)_test.json",
        )
        data_dict = LOpt.parse_file(file_path)
        augment_budget = 3
        return data_dict, augment_budget
    end

    data_dict, augment_budget = data_II()

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer)

    @test result["termination_status"] == MOI.OPTIMAL
    @test result["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result["objective"], 34.404789725, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result["solution"]["z_var"][2, 4], 1.0)
    @test isapprox(result["solution"]["z_var"][2, 5], 1.0)
end

@testset "Test hamiltonian cycle problem" begin
    function data_II()
        data_dict = Dict{String,Any}()
        data_dict["num_nodes"] = 5
        data_dict["adjacency_base_graph"] = zeros(5, 5)
        data_dict["adjacency_augment_graph"] =
            [0 10 14 12 15; 10 0 9 6 7; 14 9 0 8 9; 12 6 8 0 3; 15 7 9 3 0]
        augment_budget = 5
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_II()

    params = Dict{String,Any}(
        "data_dict" => data_dict,
        "augment_budget" => augment_budget,
        "graph_type" => "hamiltonian_cycle",
    )

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [2, data_dict["num_nodes"]],
        :projected_eigen_cuts => false,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["termination_status"] == MOI.OPTIMAL
    @test result["primal_status"] == MOI.FEASIBLE_POINT
    @test isapprox(result["objective"], 11.53910488161, atol = 1E-6)
    for (i, j) in [(1, 3), (1, 5), (2, 4), (2, 5), (3, 4)]
        @test isapprox(result["solution"]["z_var"][i, j], 1.0)
    end
end

@testset "Test topology flow/cycle and 3-minor cuts" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [3],
        :topology_flow_cuts => false,
        :projected_eigen_cuts => true,
        :solution_type => "optimal",
        :lazycuts_logging => true,
        :time_limit => test_time_limit(),
    )

    result_1 = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result_1["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_1["objective"], 11.16208322, atol = 1E-6)

    # Test flow cuts
    model_options[:max_cycle_length] = 0
    model_options[:topology_flow_cuts] = true
    result_2 = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result_2["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_2["objective"], 11.16208322, atol = 1E-6)

    # Test cycle elimination cuts
    model_options[:max_cycle_length] = num_nodes
    result_3 = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result_3["termination_status"] == MOI.OPTIMAL
    @test isapprox(result_3["objective"], 11.16208322, atol = 1E-6)
end

@testset "Test Cheeger cuts and best_incumbent warm start" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [num_nodes],
        :projected_eigen_cuts => true,
        :topology_flow_cuts => true,
        :cheeger_cuts => true,
        :best_lower_bound => 11.1591873084,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result["termination_status"] == MOI.OPTIMAL
    @test isapprox(result["objective"], 11.1591873084, atol = 1E-6)

    model_options = Dict{Symbol,Any}(:best_incumbent => [(1, 4), (2, 4), (3, 4), (3, 5)])
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result["termination_status"] == MOI.OPTIMAL
    @test isapprox(result["objective"], 11.1591873084, atol = 1E-6)
end

@testset "Test SDP relaxation" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :eigen_cuts_sizes => [2, num_nodes],
        :projected_eigen_cuts => false,
        :topology_flow_cuts => true,
        :sdp_relaxation => true,
        :solution_type => "optimal",
        :time_limit => test_time_limit(),
    )

    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result["termination_status"] == MOI.OPTIMAL
    @test isapprox(result["objective"], 23.9350775692, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1, 4], 1.0, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][1, 2], 0.425365873794, atol = 1E-6)
    @test isapprox(result["solution"]["z_var"][3, 5], 0.614568536054, atol = 1E-6)
end
