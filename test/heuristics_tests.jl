@testset "Heuristic test: k-opt" begin
    num_nodes = 5
    instance = 1
    data_dict, augment_budget = data_spanning_tree(num_nodes, instance)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => augment_budget)

    model_options =
        Dict{Symbol,Any}(:solution_type => "heuristic", :time_limit => test_time_limit())
    result = LaplacianOpt.run_LOpt(params, glpk_optimizer; options = model_options)
    @test result["solution_type"] == "heuristic"
end
