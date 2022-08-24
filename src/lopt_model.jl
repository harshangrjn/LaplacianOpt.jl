function build_LOModel(data::Dict{String,Any})
    lom = LaplacianOptModel(data)
    LOpt.variable_LOModel(lom)

    if lom.data["solution_type"] == "exact"
        LOpt.constraint_LOModel(lom)
    end

    LOpt.objective_LOModel(lom)

    return lom
end

function variable_LOModel(lom::LaplacianOptModel)
    LOpt.variable_lifted_W_matrix(lom)
    LOpt.variable_edge_onoff(lom)
    lom.data["topology_multi_commodity"] && LOpt.variable_multi_commodity_flow(lom)
    LOpt.variable_algebraic_connectivity(lom)

    return
end

function constraint_LOModel(lom::LaplacianOptModel)
    LOpt.constraint_build_W_var_matrix(lom)
    LOpt.constraint_single_vertex_cutset(lom)
    LOpt.constraint_augment_edges_budget(lom)
    lom.data["topology_multi_commodity"] &&
        LOpt.constraint_topology_multi_commodity_flow(lom)
    LOpt.constraint_lazycallback_wrapper(lom)

    return
end

function objective_LOModel(lom::LaplacianOptModel)
    LOpt.objective_maximize_algebraic_connectivity(lom)

    return
end

function optimize_LOModel!(lom::LaplacianOptModel; optimizer = nothing)
    if lom.data["relax_integrality"]
        JuMP.relax_integrality(lom.model)
    end

    JuMP.set_time_limit_sec(lom.model, lom.data["time_limit"])

    if JuMP.mode(lom.model) != JuMP.DIRECT && optimizer !== nothing
        if lom.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
            JuMP.set_optimizer(lom.model, optimizer)
        else
            Memento.warn(
                _LOGGER,
                "Model already contains optimizer, cannot use optimizer specified in `optimize_LOModel!`",
            )
        end
    end

    if JuMP.mode(lom.model) != JuMP.DIRECT &&
       lom.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
        Memento.error(
            _LOGGER,
            "No optimizer specified in `optimize_LOModel!` or the given JuMP model.",
        )
    end

    start_time = time()
    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(lom.model)

    try
        solve_time = JuMP.solve_time(lom.model)
    catch
        Memento.warn(
            _LOGGER,
            "The given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.",
        )
    end

    Memento.debug(_LOGGER, "JuMP model optimize time: $(time() - start_time)")
    lom.result = LOpt.build_LOModel_result(lom, solve_time)

    return lom.result
end

function run_LOpt(
    params::Dict{String,Any},
    lom_optimizer::MOI.OptimizerWithAttributes;
    visualize_solution = false,
    visualizing_tool = "tikz",
    display_edge_weights = false,
)
    data = LOpt.get_data(params)
    model_lopt = LOpt.build_LOModel(data)
    result_lopt = LOpt.optimize_LOModel!(model_lopt, optimizer = lom_optimizer)

    if visualize_solution
        LOpt.visualize_solution(
            result_lopt,
            data,
            visualizing_tool = visualizing_tool,
            display_edge_weights = display_edge_weights,
        )
    end

    return result_lopt
end

function run_MaxSpanTree(
    params::Dict{String,Any},
    lom_optimizer::MOI.OptimizerWithAttributes;
    visualize_solution = false,
    visualizing_tool = "tikz",
    display_edge_weights = false,
    lazy_callback = false,
)
    data = LOpt.get_data(params)
    model_mst = LOpt.build_MaxSpanTree_model(data, lazy_callback)
    result_mst = LOpt.optimize_LOModel!(model_mst, optimizer = lom_optimizer)

    if visualize_solution
        LOpt.visualize_solution(
            result_mst,
            data,
            visualizing_tool = visualizing_tool,
            display_edge_weights = display_edge_weights,
        )
    end

    return result_mst
end

function build_MaxSpanTree_model(data::Dict{String,Any}, lazy_callback::Bool)
    m_mst = LaplacianOptModel(data)
    LOpt.variable_MaxSpanTree_model(m_mst, lazy_callback)
    LOpt.constraint_MaxSpanTree_model(m_mst, lazy_callback)
    LOpt.objective_MaxSpanTree_model(m_mst)

    return m_mst
end

function objective_MaxSpanTree_model(lom::LaplacianOptModel)
    LOpt.objective_maximize_spanning_tree_cost(lom)

    return
end

function variable_MaxSpanTree_model(lom::LaplacianOptModel, lazy_callback::Bool)
    LOpt.variable_edge_onoff(lom)

    if !lazy_callback
        LOpt.variable_multi_commodity_flow(lom)
    end

    return
end

function constraint_MaxSpanTree_model(lom::LaplacianOptModel, lazy_callback::Bool)
    LOpt.constraint_single_vertex_cutset(lom)
    LOpt.constraint_augment_edges_budget(lom)

    if lazy_callback
        LOpt.constraint_lazycallback_wrapper(lom, max_span_tree = true)
    else
        LOpt.constraint_topology_multi_commodity_flow(lom)
    end

    return
end
