#------------------------------------------------------#
# Build all constraints of the LaplacianOptModel here  #
#------------------------------------------------------#

function constraint_build_W_var_matrix(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    num_edges_existing = lom.data["num_edges_existing"]
    adjacency_base_graph = lom.data["adjacency_base_graph"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    adjacency_full_graph = adjacency_augment_graph
    (num_edges_existing > 0) && (adjacency_full_graph += adjacency_base_graph)

    # Diagonal entries
    JuMP.@constraint(
        lom.model,
        [i = 1:num_nodes],
        lom.variables[:W_var][i, i] ==
        sum(adjacency_full_graph[i, :] .* lom.variables[:z_var][i, :]) -
        lom.variables[:γ_var] * (num_nodes - 1) / num_nodes
    )

    # Off-diagonal entries
    JuMP.@constraint(
        lom.model,
        [i = 1:(num_nodes-1), j = (i+1):num_nodes],
        lom.variables[:W_var][i, j] ==
        -(adjacency_full_graph[i, j] * lom.variables[:z_var][i, j]) +
        lom.variables[:γ_var] / num_nodes
    )

    return
end

function constraint_single_vertex_cutset(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    num_edges_existing = lom.data["num_edges_existing"]
    adjacency_base_graph = lom.data["adjacency_base_graph"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]
    graph_type = lom.data["graph_type"]

    adjacency_full_graph = adjacency_augment_graph
    (num_edges_existing > 0) && (adjacency_full_graph += adjacency_base_graph)

    if graph_type == "hamiltonian_cycle"
        for i in 1:num_nodes
            if sum(adjacency_full_graph[i, :] .> 0) >= 2
                JuMP.@constraint(
                    lom.model,
                    sum(lom.variables[:z_var][i, j] for j in setdiff((1:num_nodes), i)) ==
                    2
                )
            end
        end
    else # for any other graph_type
        for i in 1:num_nodes
            if sum(adjacency_full_graph[i, :] .> 0) >= 1
                JuMP.@constraint(
                    lom.model,
                    sum(lom.variables[:z_var][i, j] for j in setdiff((1:num_nodes), i)) >=
                    1
                )
            end
        end
    end

    return
end

function constraint_augment_edges_budget(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    augment_budget = lom.data["augment_budget"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    JuMP.@constraint(
        lom.model,
        sum(
            lom.variables[:z_var][i, j] for i in 1:(num_nodes-1), j in (i+1):num_nodes if
            adjacency_augment_graph[i, j] > 0
        ) == augment_budget
    )

    return
end

function constraint_lazycallback_wrapper(lom::LaplacianOptModel; optimizer = nothing)
    function constraint_lazycallback_cuts(cb_cuts)
        status = JuMP.callback_node_status(cb_cuts, lom.model)

        if status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            Memento.error(
                _LOGGER,
                "Callback status: unknown - solution may not be integer feasible",
            )

        elseif status in [MOI.CALLBACK_NODE_STATUS_INTEGER]
            #    [MOI.CALLBACK_NODE_STATUS_INTEGER, MOI.CALLBACK_NODE_STATUS_FRACTIONAL]
            if !(lom.options.formulation_type == "max_span_tree")
                if (
                    size(lom.options.eigen_cuts_sizes)[1] > 0 &&
                    minimum(lom.options.eigen_cuts_sizes) >= 2
                ) || lom.options.sdp_relaxation
                    W_val = JuMP.callback_value.(Ref(cb_cuts), lom.variables[:W_var])
                    LOpt.constraint_eigen_cuts(W_val, cb_cuts, lom)
                end

                if lom.options.soc_linearized_cuts
                    W_val = JuMP.callback_value.(Ref(cb_cuts), lom.variables[:W_var])
                    LOpt.constraint_soc_cuts_on_2minors(W_val, cb_cuts, lom)
                end

                if lom.options.cheeger_cuts
                    z_val =
                        abs.(JuMP.callback_value.(Ref(cb_cuts), lom.variables[:z_var]))
                    LOpt.round_zeros_ones!(z_val)
                    LOpt.constraint_cheeger_cuts(z_val, cb_cuts, lom, optimizer)
                end
            end

            if !(lom.data["is_base_graph_connected"]) && (lom.options.topology_flow_cuts)
                z_val = abs.(JuMP.callback_value.(Ref(cb_cuts), lom.variables[:z_var]))
                LOpt.round_zeros_ones!(z_val)
                LOpt.constraint_topology_flow_cuts(z_val, cb_cuts, lom)
            end
        end
    end

    return JuMP.set_attribute(
        lom.model,
        MOI.LazyConstraintCallback(),
        constraint_lazycallback_cuts,
    )
end

function constraint_eigen_cuts(W_val::Matrix{<:Number}, cb_cuts, lom::LaplacianOptModel)
    if (length(keys(lom.minor_idx_dict)) > 0)
        for cut_size in sort(collect(keys(lom.minor_idx_dict)), rev = true)
            for i in lom.minor_idx_dict[cut_size]
                LOpt._add_eigen_cut_lazy(W_val, cb_cuts, lom, collect(i))
            end
        end
    end
    return
end

function _add_eigen_cut_lazy(
    W_val::Matrix{<:Number},
    cb_cuts,
    lom::LaplacianOptModel,
    idx::Vector{Int64},
)
    num_nodes = lom.data["num_nodes"]

    if length(idx) in collect(2:num_nodes)
        z_var = lom.variables[:z_var]
        γ_var = lom.variables[:γ_var]

        adjacency_full_graph = lom.data["adjacency_augment_graph"]
        (lom.data["num_edges_existing"] > 0) &&
            (adjacency_full_graph += lom.data["adjacency_base_graph"])
        adjacency = adjacency_full_graph .* z_var
    end

    violated_eigen_vec =
        LOpt._violated_eigen_vector(W_val[idx, idx], tol = lom.options.tol_psd)
    if violated_eigen_vec !== nothing
        # Projected eigen-cuts in (x, gamma) space
        if lom.options.projected_eigen_cuts
            # NxN matrix
            if length(idx) == num_nodes
                con = JuMP.@build_constraint(
                    sum(
                        adjacency[i, j] *
                        (violated_eigen_vec[i] - violated_eigen_vec[j])^2 for
                        i in 1:(num_nodes-1), j in ((i+1):num_nodes)
                    ) >= γ_var
                )
            else
                # Smaller principal sub-matrices
                deg_v = Vector{JuMP.AffExpr}()
                for i in 1:num_nodes
                    push!(deg_v, sum(adjacency[i, :]))
                end

                con = JuMP.@build_constraint(
                    sum(deg_v[idx[i]] * violated_eigen_vec[i]^2 for i in 1:length(idx)) - (
                        2 * sum(
                            adjacency[idx[i], idx[j]] *
                            violated_eigen_vec[i] *
                            violated_eigen_vec[j] for i in 1:(length(idx)-1),
                            j in (i+1):length(idx)
                        )
                    ) >=
                    (γ_var / num_nodes) * (
                        num_nodes - 1 -
                        2 * sum(
                            violated_eigen_vec[i] * violated_eigen_vec[j] for
                            i in 1:(length(idx)-1), j in (i+1):length(idx)
                        )
                    )
                )
            end

        else
            # Eigen-cuts in W-space
            con = JuMP.@build_constraint(
                violated_eigen_vec' *
                lom.variables[:W_var][idx, idx] *
                violated_eigen_vec >= 0
            )
        end

        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)

        if lom.options.lazycuts_logging
            Memento.info(
                _LOGGER,
                "Polyhedral relaxation cuts: $(length(idx))x$(length(idx)) minor (eigen)",
            )
        end
    end
end

function constraint_soc_cuts_on_2minors(
    W_val::Matrix{<:Number},
    cb_cuts,
    lom::LaplacianOptModel,
)
    num_nodes = lom.data["num_nodes"]

    for i in 1:(num_nodes-1)
        for j in (i+1):num_nodes
            W_val_minor = W_val[i, j]^2 - W_val[i, i] * W_val[j, j]
            if W_val_minor >= lom.options.tol_zero
                coefficient_vec = [
                    -(W_val[i, j])^2
                    2 * W_val[i, j] * W_val[i, i]
                ]

                variable_vec = [lom.variables[:W_var][i, i], lom.variables[:W_var][i, j]]

                # Taylor's approximation of rotated SOC constraint (x^2 <= y*z)
                if !isapprox(W_val[i, i], 0, atol = 1E-6)
                    con = JuMP.@build_constraint(
                        coefficient_vec[1] * variable_vec[1] +
                        coefficient_vec[2] * variable_vec[2] <=
                        W_val[i, i]^2 * lom.variables[:W_var][j, j]
                    )

                    MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)
                end

                if lom.options.lazycuts_logging
                    Memento.info(
                        _LOGGER,
                        "Polyhedral relaxation cuts: 2x2 minor linearized SOC",
                    )
                end
            end
        end
    end

    return
end

function constraint_cheeger_cuts(
    z_val::Matrix{<:Number},
    cb_cuts,
    lom::LaplacianOptModel,
    optimizer::MOI.OptimizerWithAttributes,
)
    result_dict =
        LOpt.cheeger_constant(z_val .* lom.data["adjacency_augment_graph"], optimizer)

    ac_lower_bound = floor(lom.options.best_lower_bound, digits = 3)

    if !(isapprox(result_dict["cheeger_constant"], 0, atol = 1E-5)) && (
        result_dict["cheeger_constant"] <
        (ac_lower_bound * lom.options.cheeger_cuts_factor)
    )
        lom.result["cheeger_cuts_separation_time"] += result_dict["solve_time"]
        S = result_dict["S"]
        S_complement = result_dict["S_complement"]

        con = JuMP.@build_constraint(
            sum(
                lom.data["adjacency_augment_graph"][i, j] * lom.variables[:z_var][i, j] for i in S, j in S_complement
            ) >=
            ac_lower_bound * lom.options.cheeger_cuts_factor * length(collect(S))
        )
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)
    end
end

function constraint_topology_flow_cuts(
    z_val::Matrix{<:Number},
    cb_cuts,
    lom::LaplacianOptModel,
)
    num_nodes = lom.data["num_nodes"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]
    graph_type = lom.data["graph_type"]
    num_edges_existing = lom.data["num_edges_existing"]
    augment_budget = lom.data["augment_budget"]

    G = Graphs.SimpleGraph(abs.(z_val))
    cc_lazy = Graphs.connected_components(G)

    max_cc = 5 # increase this to any greater integer value and it works
    if length(cc_lazy) > max_cc
        Memento.info(
            _LOGGER,
            "Polyhedral relaxation: flow cuts not added for integer solutions with $(length(cc_lazy)) connected components",
        )
    elseif length(cc_lazy) == 1
        return
    end

    is_spanning_tree = false
    if (num_edges_existing == 0) && (augment_budget == (num_nodes - 1))
        is_spanning_tree = true
    end

    # Subtour elimination (for cycle/tree graphs)
    cycle_elimination = false

    if ((graph_type == "hamiltonian_cycle") || (is_spanning_tree))
        cycles = Set()
        (lom.options.max_cycle_length >= 3) &&
            (cycles = LOpt.get_unique_cycles(G, lom.options.max_cycle_length))

        if length(cycles) >= 1
            cycle_elimination = true

            for i in 1:length(cycles)
                cyc = collect(cycles)[i]
                cyc_length = length(cyc)
                expr = JuMP.AffExpr()

                for ii in 1:cyc_length
                    jj = ii + 1
                    (ii == cyc_length) && (jj = 1)
                    if !(isapprox(
                        adjacency_augment_graph[cyc[ii], cyc[jj]],
                        0,
                        atol = 1E-6,
                    ))
                        expr += lom.variables[:z_var][
                            min(cyc[ii], cyc[jj]),
                            max(cyc[ii], cyc[jj]),
                        ]
                    end
                end
                con = JuMP.@build_constraint(expr <= (cyc_length - 1))
                MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)

                if lom.options.lazycuts_logging
                    Memento.info(
                        _LOGGER,
                        "Polyhedral relaxation cuts: cycle elimination added (size = $cyc_length)",
                    )
                end
            end
        end
    end

    # Flow cuts on connected components    
    if !(cycle_elimination) && (length(cc_lazy) in 2:max_cc)
        num_edges_cutset = 1
        if graph_type == "hamiltonian_cycle"
            num_edges_cutset = 2
        end

        for k in 1:(length(cc_lazy)-1)
            # From cutset  
            cutset_f = cc_lazy[k]
            # To cutset
            ii = setdiff(1:length(cc_lazy), k)
            cutset_t = reduce(vcat, cc_lazy[ii])
            if LOpt._is_flow_cut_valid(
                cutset_f,
                cutset_t,
                num_edges_cutset,
                adjacency_augment_graph,
            )
                con = JuMP.@build_constraint(
                    sum(
                        lom.variables[:z_var][i, j] for
                        i in cutset_f, j in cutset_t if
                        !(isapprox(adjacency_augment_graph[i, j], 0, atol = 1E-6))
                    ) >= num_edges_cutset
                )

                MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)
            end

            if lom.options.lazycuts_logging
                Memento.info(_LOGGER, "Polyhedral relaxation cuts: flow cuts added")
            end
        end
    end

    return
end

function constraint_topology_multi_commodity_flow(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]

    # Flow balance on source node 
    JuMP.@constraint(
        lom.model,
        [k = 1:(num_nodes-1)],
        sum(lom.variables[:flow_var][1, 2:num_nodes, k]) -
        sum(lom.variables[:flow_var][2:num_nodes, 1, k]) == 1
    )
    # Flow balance on target node 
    JuMP.@constraint(
        lom.model,
        [k = 1:(num_nodes-1)],
        sum(lom.variables[:flow_var][:, k+1, k]) -
        sum(lom.variables[:flow_var][k+1, :, k]) == 1
    ) #(k-th commodity corresponds to k+1th node)
    # Flow balance on remaining nodes 
    JuMP.@constraint(
        lom.model,
        [k = 1:(num_nodes-1), v = (2:num_nodes); k != v - 1],
        sum(lom.variables[:flow_var][:, v, k]) - sum(lom.variables[:flow_var][v, :, k]) ==
        0
    )
    # Flow <= capacity 
    JuMP.@constraint(
        lom.model,
        [k = 1:(num_nodes-1), i = 1:num_nodes, j = 1:num_nodes, kdash = 1:(num_nodes-1)],
        lom.variables[:flow_var][i, j, k] + lom.variables[:flow_var][j, i, kdash]' <=
        lom.variables[:z_var][i, j]
    )

    return
end

function constraint_sdp_relaxation_dummy(lom::LaplacianOptModel)
    JuMP.@constraint(lom.model, sum(lom.variables[:sdp_dummy_var]) == 1.0)

    return
end
