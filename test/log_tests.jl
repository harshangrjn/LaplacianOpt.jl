@testset "visualize_solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/optimizer.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "tol_zero" => 1E-4,
        "tol_psd" => 1E-3,
        "eigen_cuts_full" => true,
        "topology_flow_cuts" => true,
        "lazycuts_logging" => true
    )

    data = LOpt.get_data(params)
    model_lopt  = LOpt.build_LOModel(data)
    result_lopt = LOpt.optimize_LOModel!(model_lopt, optimizer = glpk_optimizer)

    LOpt.visualize_solution(result_lopt, data, visualizing_tool = "tikz", plot_file_format = "tex")
    LOpt.visualize_solution(result_lopt, data, visualizing_tool = "tikz", plot_file_format = "tex", display_edge_weights=true)
    LOpt.visualize_solution(result_lopt, data, visualizing_tool = "graphviz")
    LOpt.visualize_solution(result_lopt, data, visualizing_tool = "graphviz", display_edge_weights=false)

    @test result_lopt["termination_status"] == MOI.OPTIMAL
    @test result_lopt["primal_status"] == MOI.FEASIBLE_POINT

    # Nothing more to test for coverage since the visualizaiton depends on exterior packages. So, I believe the dot files generated are accurate. 
end