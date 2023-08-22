function build_LOModel(data::Dict{String,Any}; optimizer = nothing, options = nothing)
    lom = LOpt.LaplacianOptModel(data)

    # Default eigen-cuts sizes (for exact soltution of MISDP, include `num_nodes` here)
    LOpt.set_option(lom, :eigen_cuts_sizes, [lom.data["num_nodes"], 2])

    # Update defaults to user-defined options
    if !isnothing(options)
        for i in keys(options)
            LOpt.set_option(lom, i, options[i])
        end
    end

    if lom.options.formulation_type == "max_λ2"
        if lom.options.solution_type == "heuristic"
            LOpt._logging_info(lom)
            return lom

        elseif lom.options.solution_type == "optimal"
            # Populate PMinor indices
            lom.minor_idx_dict =
                LOpt._PMinorIdx(lom.data["num_nodes"], lom.options.eigen_cuts_sizes)
            LOpt._logging_info(lom)

            LOpt.variable_LOModel(lom)
            LOpt.constraint_LOModel(lom; optimizer = optimizer)
            LOpt.objective_LOModel(lom)
        end
    elseif lom.options.formulation_type == "max_span_tree"
        if lom.options.solution_type in ["optimal"]
            LOpt.variable_MaxSpanTree_model(lom)
            LOpt.constraint_MaxSpanTree_model(lom)
            LOpt.objective_MaxSpanTree_model(lom)
        end
    end

    return lom
end

function variable_LOModel(lom::LaplacianOptModel)
    LOpt.variable_lifted_W_matrix(lom)
    LOpt.variable_edge_onoff(lom)
    if !lom.data["is_base_graph_connected"] && lom.options.topology_multi_commodity
        LOpt.variable_multi_commodity_flow(lom)
    end
    LOpt.variable_algebraic_connectivity(lom)
    lom.options.sdp_relaxation && LOpt.variable_sdp_relaxation_dummy(lom)

    return
end

function constraint_LOModel(lom::LaplacianOptModel; optimizer = nothing)
    LOpt.constraint_build_W_var_matrix(lom)
    LOpt.constraint_single_vertex_cutset(lom)
    LOpt.constraint_augment_edges_budget(lom)
    if !lom.data["is_base_graph_connected"] && lom.options.topology_multi_commodity
        LOpt.constraint_topology_multi_commodity_flow(lom)
    end
    lom.options.sdp_relaxation && LOpt.constraint_sdp_relaxation_dummy(lom)
    LOpt.lazycallback_status(lom) &&
        LOpt.constraint_lazycallback_wrapper(lom, optimizer = optimizer)

    return
end

function objective_LOModel(lom::LaplacianOptModel)
    LOpt.objective_maximize_algebraic_connectivity(lom)

    return
end

function optimize_LOModel!(lom::LaplacianOptModel; optimizer = nothing)
    if lom.options.relax_integrality
        JuMP.relax_integrality(lom.model)
    end

    JuMP.set_time_limit_sec(lom.model, lom.options.time_limit)

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
    result_dict = LOpt.build_LOModel_result(lom, solve_time)

    return merge(lom.result, result_dict)
end

function run_LOpt(
    params::Dict{String,Any},
    lom_optimizer::MOI.OptimizerWithAttributes;
    visualize_solution = false,
    visualizing_tool = "tikz",
    display_edge_weights = false,
    options = nothing,
)
    data = LOpt.get_data(params)
    model_lopt = LOpt.build_LOModel(data, optimizer = lom_optimizer, options = options)

    if model_lopt.options.solution_type == "heuristic"
        result_lopt = LOpt.heuristic_kopt(model_lopt)
    elseif model_lopt.options.solution_type == "optimal"
        result_lopt = LOpt.optimize_LOModel!(model_lopt, optimizer = lom_optimizer)
    end

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

function objective_MaxSpanTree_model(lom::LaplacianOptModel)
    LOpt.objective_maximize_spanning_tree_cost(lom)

    return
end

function variable_MaxSpanTree_model(lom::LaplacianOptModel)
    LOpt.variable_edge_onoff(lom)

    if lom.options.topology_multi_commodity
        LOpt.variable_multi_commodity_flow(lom)
    end

    return
end

function constraint_MaxSpanTree_model(lom::LaplacianOptModel)
    LOpt.constraint_single_vertex_cutset(lom)
    LOpt.constraint_augment_edges_budget(lom)

    if lom.options.topology_flow_cuts
        LOpt.constraint_lazycallback_wrapper(lom)
    elseif lom.options.topology_multi_commodity
        LOpt.constraint_topology_multi_commodity_flow(lom)
    end

    return
end

function set_option(lom::LaplacianOptModel, s::Symbol, val)
    return Base.setproperty!(lom.options, s, val)
end

function lazycallback_status(lom::LaplacianOptModel)
    if (
           size(lom.options.eigen_cuts_sizes)[1] > 0 &&
           minimum(lom.options.eigen_cuts_sizes) >= 2
       ) ||
       lom.options.topology_flow_cuts ||
       lom.options.soc_linearized_cuts ||
       lom.options.cheeger_cuts ||
       lom.options.sdp_relaxation
        return true
    else
        return false
    end
end

function _logging_info(lom::LaplacianOptModel)
    if lom.options.solution_type == "optimal"
        if length(keys(lom.minor_idx_dict)) > 0
            for k in keys(lom.minor_idx_dict)
                Memento.info(_LOGGER, "Applying eigen cuts ($(k)x$(k) matrix)")
            end
        end

        lom.options.soc_linearized_cuts &&
            Memento.info(_LOGGER, "Applying linearized SOC cuts (2x2 minors)")
        lom.options.cheeger_cuts &&
            Memento.info(_LOGGER, "Applying Cheeger inequality-based cuts")
        lom.options.sdp_relaxation &&
            Memento.info(_LOGGER, "Applying SDP relaxation formulation")

        if !lom.data["is_base_graph_connected"] && lom.options.topology_multi_commodity
            Memento.info(_LOGGER, "Applying topology multi commodity constraints")
        end

        if !lom.data["is_base_graph_connected"] && lom.options.topology_flow_cuts
            Memento.info(_LOGGER, "Applying topology flow cuts")
        end
    elseif lom.options.solution_type == "heuristic"
        Memento.info(_LOGGER, "Applying heuristics to obtain a feasible solution")
    end
end
