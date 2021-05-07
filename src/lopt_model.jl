function build_LOModel(data::Dict{String, Any})

    m_lo = LaplacianOptModel(data, JuMP.Model(), Dict{Symbol,Any}(), Dict{String,Any}())

    variable_LOModel(m_lo)
    
    if m_lo.data["solution_type"] == "exact"
        constraint_LOModel(m_lo)
    end

    objective_LOModel(m_lo)

    return m_lo
end

function variable_LOModel(lom::LaplacianOptModel)

    variable_lifted_W_matrix(lom)
    variable_edge_onoff(lom)
    variable_algebraic_connectivity(lom)

    return
end

function constraint_LOModel(lom::LaplacianOptModel)
    
    constraint_build_W_var_matrix(lom)

    constraint_topology_no_self_loops(lom)
    constraint_topology_vertex_cutset(lom)
    constraint_topology_total_edges(lom)

    constraint_lazycallback_wrapper(lom)
    
    return
end

function objective_LOModel(lom::LaplacianOptModel)
    objective_maximize_algebraic_connectivity(lom)

    return
end

function optimize_LOModel!(lom::LaplacianOptModel; optimizer=nothing)
    if lom.data["relax_integrality"]
        JuMP.relax_integrality(lom.model)
    end

    if JuMP.mode(lom.model) != JuMP.DIRECT && optimizer !== nothing
        if lom.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
            JuMP.set_optimizer(lom.model, optimizer)
        else
            Memento.warn(_LOGGER, "Model already contains optimizer, cannot use optimizer specified in `optimize_LOModel!`")
        end
    end

    if JuMP.mode(lom.model) != JuMP.DIRECT && lom.model.moi_backend.state == MOI.Utilities.NO_OPTIMIZER
        Memento.error(_LOGGER, "No optimizer specified in `optimize_LOModel!` or the given JuMP model.")
    end
    
    start_time = time()

    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(lom.model)

    try
        solve_time = JuMP.solve_time(lom.model)
    catch
        Memento.warn(_LOGGER, "The given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.");
    end
    
    Memento.debug(_LOGGER, "JuMP model optimize time: $(time() - start_time)")
    
    lom.result = build_LOModel_result(lom, solve_time) 

    return lom.result
end

function run_LOpt_model(params::Dict{String, Any}, lom_optimizer::MOI.OptimizerWithAttributes; visualize_solution = false, visualizing_tool = "graphviz")

    data = LO.get_data(params)

    model_lopt  = LO.build_LOModel(data)

    result_lopt = LO.optimize_LOModel!(model_lopt, optimizer = lom_optimizer)

    if visualize_solution
        LaplacianOpt.visualize_LOModel_solution(result_lopt, data, visualizing_tool = visualizing_tool)
    end

    return result_lopt
end