@testset "visualize_LOModel_solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/solver.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "optimizer" => "glpk",
        "tol_zero" => 1E-4,
        "tol_psd" => 1E-3,
        "eigen_cuts_full" => true,
        "topology_flow_cuts" => true,
        "lazycuts_logging" => true
    )

    lom_optimizer = get_solver(params)
    data = LO.get_data(params)
    model_lopt  = LO.build_LOModel(data)
    result_lopt = LO.optimize_LOModel!(model_lopt, optimizer = lom_optimizer)

    LO.visualize_LOModel_solution(result_lopt, data, visualizing_tool = "tikz", plot_file_format = "tex")
    LO.visualize_LOModel_solution(result_lopt, data, visualizing_tool = "graphviz")

    @test result_lopt["termination_status"] == MOI.OPTIMAL
    @test result_lopt["primal_status"] == MOI.FEASIBLE_POINT

    # Nothing more to test for coverage since the visualizaiton depends on exterior packages. So, we believe the dot files generated are accurate. 
end