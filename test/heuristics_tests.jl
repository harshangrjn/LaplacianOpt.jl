@testset "Heuristic test: 1-opt_spanning_tree" begin
    num_nodes = 8
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
            :solution_type => "heuristic", 
            :kopt_parameter => 1,
            :num_central_nodes_verifier => 5,
            :time_limit => test_time_limit()
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 22.8042, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 6], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][6, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 8], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][8, 7], 1.0)
end

@testset "Heuristic test: 2-opt_spanning_tree" begin
    num_nodes = 8
    instance = 5
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 2,
        :num_central_nodes_verifier => 5,
        :time_limit => test_time_limit(),
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 22.5051, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 6], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][6, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 7], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][7, 8], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][8, 7], 1.0)

end

@testset "Heuristic test: 3-opt_spanning_tree" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 3,
        :num_central_nodes_verifier => 5,
        :time_limit => test_time_limit(),
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 11.1592, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 3], 1.0)
end

@testset "Heuristic test: 1-opt_tree" begin
    num_nodes = 5
    instance = 11
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 1,
        :num_central_nodes_verifier => 3,
        :time_limit => test_time_limit(),
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 36.0425, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 4], 1.0)
end

@testset "Heuristic test: 2-opt_tree" begin
    num_nodes = 5
    instance = 11
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 2,
        :num_central_nodes_verifier => 3,
        :time_limit => test_time_limit(),
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 36.0425, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 4], 1.0)
end

@testset "Heuristic test: 3-opt_tree" begin
    num_nodes = 5
    instance = 11
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 3,
        :num_central_nodes_verifier => 3,
        :time_limit => test_time_limit(),
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result_lo["objective"], 36.0425, atol = 1E-4)
    @test isapprox(result_lo["solution"]["z_var"][1, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][1, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 1], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][2, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 2], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 4], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][3, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 3], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][4, 5], 1.0)
    @test isapprox(result_lo["solution"]["z_var"][5, 4], 1.0)
end
