@testset "visualize_LOModel_solution tests" begin
    
    file_path = joinpath(@__DIR__,"..", "examples/solver.jl")
    include(file_path)

    params = Dict{String, Any}(
        "num_nodes" => 5,
        "instance" => 1,
        "optimizer" => "glpk"
    )

    lom_optimizer = get_solver(params)
    data = LO.get_data(params)
    model_lopt  = LO.build_LOModel(data)
    result_lopt = LO.optimize_LOModel!(model_lopt, optimizer = lom_optimizer)

    LO.visualize_LOModel_solution(result_lopt, data, visualizing_tool = "tikz")
    LO.visualize_LOModel_solution(result_lopt, data, visualizing_tool = "graphviz")
    
end