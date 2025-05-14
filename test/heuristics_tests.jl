@testset "Heuristic test: 1-opt_spanning_tree" begin
    num_nodes = 8
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 1,
        :num_central_nodes_kopt => 5,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 22.8042, atol = 1E-4)
    
    for (i, j) in [(1,7), (2,7), (3,7), (4,6), (4,7), (5,7), (7,8)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
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
        :num_central_nodes_kopt => 5,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 22.5051, atol = 1E-4)
    
    for (i, j) in [(1,5), (2,5), (3,5), (4,5), (5,6), (5,7), (7,8)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
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
        :num_central_nodes_kopt => 5,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 11.1592, atol = 1E-4)
    
    for (i, j) in [(1,4), (2,4), (3,4), (3,5)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
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
        :num_central_nodes_kopt => 3,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 46.028980895, atol = 1E-4)
    
    for (i, j) in [(1,2), (1,3), (1,4), (2,4), (2,5), (3,4), (3,5), (4,5)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
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
        :num_central_nodes_kopt => 3,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 46.028980895, atol = 1E-4)
    
    for (i, j) in [(1,2), (1,3), (1,4), (2,4), (2,5), (3,4), (3,5), (4,5)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
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
        :num_central_nodes_kopt => 3,
    )
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 46.028980895, atol = 1E-4)
    
    for (i, j) in [(1,2), (1,3), (1,4), (2,4), (2,5), (3,4), (3,5), (4,5)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
end

@testset "Heuristic test: zero augment budget" begin
    function data_II()
        data_dict = Dict{String,Any}()
        data_dict["num_nodes"] = 4
        data_dict["adjacency_base_graph"] = [0 2 0 0; 2 0 3 0; 0 3 0 4; 0 0 4 0]
        data_dict["adjacency_augment_graph"] = [0 0 4 8; 0 0 0 7; 4 0 0 0; 8 7 0 0]
        augment_budget = 0
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_II()

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(:solution_type => "heuristic", :kopt_parameter => 1)
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 1.610684749, atol = 1E-4)
    # Check for expected edges in the solution
    for (i, j) in [(1,2), (2,3), (3,4)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
end

@testset "Heuristic test: refinement_tree_1opt_2opt!" begin
    function data_II()
        data_dict = Dict{String,Any}()
        data_dict["num_nodes"] = 4
        data_dict["adjacency_base_graph"] = [0 2 0 0; 2 0 3 0; 0 3 0 4; 0 0 4 0]
        data_dict["adjacency_augment_graph"] = [0 0 4 8; 0 0 0 7; 4 0 0 0; 8 7 0 0]
        augment_budget = 2
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_II()

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(:solution_type => "heuristic", :kopt_parameter => 2)
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 7.8551986404, atol = 1E-4)
    
    for (i, j) in [(1,2), (2,3), (3,4), (1,4), (2,4)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
end

@testset "Heuristic test: refinement_tree_2opt_3opt!" begin
    function data_II()
        data_dict = Dict{String,Any}()
        data_dict["num_nodes"] = 5
        data_dict["adjacency_base_graph"] =
            [0 2 0 0 0; 2 0 3 0 0; 0 3 0 4 0; 0 0 4 0 5; 0 0 0 5 0]
        data_dict["adjacency_augment_graph"] =
            [0 0 4 8 10; 0 0 0 7 9; 4 0 0 0 6; 8 7 0 0 0; 10 9 6 0 0]
        augment_budget = 3
        return data_dict, augment_budget
    end
    data_dict, augment_budget = data_II()

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options = Dict{Symbol,Any}(:solution_type => "heuristic", :kopt_parameter => 3)
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)

    @test result["solution_type"] == "heuristic"
    @test isapprox(result["heuristic_objective"], 9.276487957, atol = 1E-4)
    
    for (i, j) in [(1,2), (2,3), (3,4), (4,5), (1,3), (1,5), (2,5)]
        @test isapprox(result["heuristic_solution"][i, j], 1.0)
        @test isapprox(result["heuristic_solution"][j, i], 1.0)
    end
end
