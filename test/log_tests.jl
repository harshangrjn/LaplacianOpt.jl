@testset "visualize_solution tests" begin
    num_nodes = 5
    instance = 1
    file_path = joinpath(
        @__DIR__,
        "..",
        "examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json",
    )
    data_dict = LOpt.parse_file(file_path)

    params =
        Dict{String,Any}("data_dict" => data_dict, "augment_budget" => (num_nodes - 1))

    model_options = Dict{Symbol,Any}(
        :tol_zero => 1E-4,
        :tol_psd => 1E-3,
        :eigen_cuts_full => true,
        :topology_flow_cuts => true,
        :lazycuts_logging => true,
        :time_limit => test_time_limit(),
    )

    data = LOpt.get_data(params)
    model_lopt = LOpt.build_LOModel(data)
    result_lopt = LOpt.optimize_LOModel!(model_lopt, optimizer = glpk_optimizer)

    LOpt.visualize_solution(
        result_lopt,
        data,
        visualizing_tool = "tikz",
        plot_file_format = "tex",
    )
    LOpt.visualize_solution(
        result_lopt,
        data,
        visualizing_tool = "tikz",
        plot_file_format = "tex",
        display_edge_weights = true,
    )
    LOpt.visualize_solution(result_lopt, data, visualizing_tool = "graphviz")
    LOpt.visualize_solution(
        result_lopt,
        data,
        visualizing_tool = "graphviz",
        display_edge_weights = false,
    )

    @test result_lopt["termination_status"] == MOI.OPTIMAL
    @test result_lopt["primal_status"] == MOI.FEASIBLE_POINT

    # Visualizing based on heuristic solution
    model_options = Dict{Symbol,Any}(
        :solution_type => "heuristic",
        :kopt_parameter => 2,
        :time_limit => test_time_limit(),
    )
    result_heuristic = LaplacianOpt.run_LOpt(
        params,
        glpk_optimizer;
        options = model_options,
        visualize_solution = true,  # Make it true to plot the graph solution
        visualizing_tool = "tikz",
    )
    @test result_heuristic["solution_type"] == "heuristic"
    @test isapprox(
        result_heuristic["heuristic_objective"],
        result_lopt["objective"],
        atol = 1E-4,
    )
end
