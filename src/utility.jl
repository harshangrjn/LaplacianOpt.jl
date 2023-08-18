"""
    round_zeros_ones!(v::Array{<:Number}; tol = 1E-6)

Given a vector of numbers, this function updates the input vector by rounding the values closest to 0, 1 and -1. 
"""
function round_zeros_ones!(v::Array{<:Number}; tol = 1E-6)
    for i in findall(abs.(v) .<= tol)
        v[i] = 0
    end

    for i in findall(((1 - tol) .<= abs.(v) .<= (1 + tol)))
        if v[i] > 0
            v[i] = 1
        elseif v[i] < 0
            v[i] = -1
        end
    end
end

"""
    optimal_graph_edges(adjacency_matrix::Array{<:Number}) 

Returns a vector of tuples of edges corresponding 
to an input adjacency matrix of the graph. 
"""
function optimal_graph_edges(adjacency_matrix::Array{<:Number})
    edges = Vector{Tuple{Int64,Int64}}()

    num_nodes = size(adjacency_matrix)[1]

    if !LA.issymmetric(adjacency_matrix)
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")
    end

    for i in 1:num_nodes
        if !isapprox(adjacency_matrix[i, i], 0)
            Memento.error(_LOGGER, "Input adjacency matrix cannot have self loops")
        end
    end

    for i in 1:(num_nodes-1)
        for j in (i+1):num_nodes
            if !isapprox(adjacency_matrix[i, j], 0)
                push!(edges, (i, j))
            end
        end
    end

    return edges
end

"""
    laplacian_matrix(adjacency_matrix::Array{<:Number}) 

Given a weighted adjacency matrix as an input, this function returns the 
weighted Laplacian matrix of the graph. 
"""
function laplacian_matrix(adjacency_matrix::Array{<:Number})
    if !LA.issymmetric(adjacency_matrix)
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")
    end

    num_nodes = size(adjacency_matrix)[1]

    laplacian_matrix = zeros(Float64, num_nodes, num_nodes)

    for i in 1:num_nodes
        for j in i:num_nodes
            if i == j
                if !isapprox(adjacency_matrix[i, j], 0)
                    Memento.error(
                        _LOGGER,
                        "Input adjacency matrix cannot have self loops",
                    )
                end

                laplacian_matrix[i, j] = sum(adjacency_matrix[i, :])
            else
                if adjacency_matrix[i, j] <= -1E-6
                    Memento.error(
                        _LOGGER,
                        "Input adjacency matrix cannot have negative weights",
                    )
                end

                if !isapprox(adjacency_matrix[i, j], 0)
                    laplacian_matrix[i, j] = -adjacency_matrix[i, j]
                    laplacian_matrix[j, i] = -adjacency_matrix[i, j]
                end
            end
        end
    end

    return laplacian_matrix
end

"""
    fiedler_vector(adjacency_matrix::Array{<:Number})

Returns the Fiedler vector or the eigenvector corresponding to the 
second smallest eigenvalue of the Laplacian matrix for an input weighted adjacency matrix of the graph. 
"""
function fiedler_vector(adjacency_matrix::Array{<:Number})
    L_mat = LOpt.laplacian_matrix(adjacency_matrix)
    ac = LOpt.algebraic_connectivity(adjacency_matrix)

    fiedler = LA.eigvecs(L_mat)[:, 2]

    if !isapprox(fiedler' * L_mat * fiedler, ac, atol = 1E-6)
        Memento.error(_LOGGER, "Evaluated eigenvector is not the Fiedler vector")
    end

    return fiedler
end

"""
    algebraic_connectivity(adjacency_matrix::Array{<:Number}) 
    
Returns the algebraic connectivity or the  
second smallest eigenvalue of the Laplacian matrix, for an input weighted adjacency matrix of the graph. 
"""
function algebraic_connectivity(adjacency_matrix::Array{<:Number})
    L_mat = LOpt.laplacian_matrix(adjacency_matrix)

    ac = sort(LA.eigvals(L_mat))[2]
    return ac
end

"""
    _is_flow_cut_valid(cutset_f::Vector{Int64}, cutset_t::Vector{Int64}, adjacency::Array{<:Number})

Given from and to sets of vertices of connected components of a graph, this function returns a 
boolean if the input `num_edges` number of candidate edge exist to augment or not, between these two sets of vertices. 
"""
function _is_flow_cut_valid(
    cutset_f::Vector{Int64},
    cutset_t::Vector{Int64},
    num_edges::Int64,
    adjacency::Array{<:Number},
)
    k = 0
    for i in cutset_f
        for j in cutset_t
            if !(isapprox(adjacency[i, j], 0, atol = 1E-6))
                k += 1
                if k == num_edges
                    return true
                end
            end
        end
    end

    return false
end

"""
    _violated_eigen_vector(W::Array{<:Number}; tol = 1E-6)
    
Given a square symmetric matrix, this function returns an eigen vector corresponding to it's most 
violated eigen value w.r.t positive semi-definiteness of the input matrix.  
"""
function _violated_eigen_vector(W::Array{<:Number}; tol = 1E-6)
    W_eigvals = LA.eigvals(W)

    if typeof(W_eigvals) != Vector{Float64}
        if !(isapprox(imag(W_eigvals), zeros(size(W_eigvals)[1]), atol = 1E-6))
            Memento.error(_LOGGER, "PSD matrix (W) cannot have complex eigenvalues")
        end
    end

    if LA.eigmin(W) <= -tol
        if W_eigvals[1] == LA.eigmin(W)
            violated_eigen_vec = LA.eigvecs(W)[:, 1]
            LOpt.round_zeros_ones!(violated_eigen_vec)
            return violated_eigen_vec
        else
            Memento.warn(
                _LOGGER,
                "Eigen cut corresponding to the negative eigenvalue could not be evaluated",
            )
        end
    else
        return
    end
end

"""
    cheeger_constant(G::Matrix{<:Number}, 
                     optimizer::MOI.OptimizerWithAttributes; 
                     cheeger_ub = 1E6,
                     optimizer_log = false)
    
Given a symmetric, weighted adjacency matrix `G` of a connected graph, and a mixed-integer linear optimizer, 
this function returns the Cheeger's constant (or isoperimetric number of the graph) and the associated two partitions (`S` and `S_complement`) of the graph.
References: 
(I) Chung, F.R., "Laplacians of graphs and Cheeger's inequalities",
Combinatorics, Paul Erdos is Eighty, 2(157-172), pp.13-2, 1996.
(II) Mohar, B., 1989. Isoperimetric numbers of graphs. Journal of combinatorial theory, 
Series B, 47(3), pp.274-291. Link: https://doi.org/10.1016/0095-8956(89)90029-4
"""
function cheeger_constant(
    G::Matrix{<:Number},
    optimizer::MOI.OptimizerWithAttributes;
    cheeger_ub = 1E4,
    optimizer_log = false,
)
    if size(G)[1] !== size(G)[2]
        Memento.error(_LOGGER, "Input adjacency has to be a square matrix")
    end
    num_nodes = size(G)[1]

    result_cheeger = Dict{String,Any}("cheeger_constant" => 0)

    if !isapprox(LOpt.algebraic_connectivity(G), 0, atol = 1E-5)
        m_cheeger = JuMP.Model(optimizer)
        if !optimizer_log
            JuMP.set_silent(m_cheeger)
        end

        # Variables
        JuMP.@variable(m_cheeger, 0 <= z[1:num_nodes] <= 1, Bin)
        JuMP.@variable(m_cheeger, -1 <= Z[1:num_nodes, 1:num_nodes] <= 1)
        JuMP.@variable(m_cheeger, 0 <= cheeger <= cheeger_ub)
        JuMP.@variable(m_cheeger, cheeger_z[1:num_nodes] >= 0)

        # Constraints
        JuMP.@constraint(m_cheeger, sum(z) >= 1)
        JuMP.@constraint(m_cheeger, sum(z) <= floor(num_nodes / 2))
        JuMP.@constraint(
            m_cheeger,
            sum(cheeger_z) >= sum(
                G[i, j] * Z[i, j] for i in 1:(num_nodes-1), j in (i+1):num_nodes if
                !(isapprox(G[i, j], 0, atol = 1E-5))
            )
        )

        JuMP.@constraint(
            m_cheeger,
            [i = 1:(num_nodes-1), j = (i+1):num_nodes],
            Z[i, j] >= z[i] - z[j]
        )
        JuMP.@constraint(
            m_cheeger,
            [i = 1:(num_nodes-1), j = (i+1):num_nodes],
            Z[i, j] >= z[j] - z[i]
        )

        for i in 1:num_nodes
            LOpt.relaxation_bilinear(m_cheeger, cheeger_z[i], cheeger, z[i])
        end

        # Objective
        JuMP.@objective(m_cheeger, Min, cheeger)

        JuMP.optimize!(m_cheeger)

        result_cheeger["cheeger_constant"] = JuMP.objective_value(m_cheeger)
        result_cheeger["S"] =
            Set(findall(isapprox.(abs.(JuMP.value.(z)), 1, atol = 1E-5)))
        result_cheeger["S_complement"] =
            Set(findall(isapprox.(abs.(JuMP.value.(z)), 0, atol = 1E-5)))
        result_cheeger["solve_time"] = JuMP.solve_time(m_cheeger)
    end

    return result_cheeger
end

"""
    get_unique_cycles(G::Graphs.SimpleGraphs.SimpleGraph{<:Number}, 
        max_cycle_length::Int, 
        min_cycle_length = 3, 
        ceiling = 100)

Given an undirected graph of type `Graphs.SimpleGraphs`, 
this function returns unique cycles of maximum length given by `max_cycle_length`. 
"""
function get_unique_cycles(
    G::Graphs.SimpleGraphs.SimpleGraph{<:Number},
    max_cycle_length::Int,
    min_cycle_length = 3,
    ceiling = 100,
)
    cycles = Graphs.simplecycles_limited_length(G, max_cycle_length, ceiling)
    cycles_set = Set()
    for i in 1:length(cycles)
        if length(cycles[i]) >= min_cycle_length
            push!(cycles_set, cycles[i])
        end
    end

    return collect(cycles_set)[unique(i -> sort.(cycles_set)[i], 1:length(cycles_set))]
end

"""
    get_minor_idx(num_nodes::Int64, size::Int64)

Given the number of nodes (size of the square matrix) and the 
preferred size of principal sub-matrices (`<= num_nodes`), this function 
outputs all the indices of correspondingly sized sub-matrices.  
"""
function get_minor_idx(num_nodes::Int64, size::Int64)
    minor_idx = Vector{Tuple{Int64,Vararg{Int64}}}()

    for i in 1:num_nodes
        push!(minor_idx, (i,))
    end

    if size >= 2
        for _ in 1:(size-1)
            global minor_idx_tmp = Vector{Tuple{Int64,Vararg{Int64}}}()
            for k in 1:length(minor_idx)
                idx_k = minor_idx[k]
                if idx_k[end] <= (num_nodes - 1)
                    for ii in collect((idx_k[end]+1):num_nodes)
                        push!(minor_idx_tmp, Tuple(push!(collect(idx_k), ii)))
                    end
                end
            end
            minor_idx = deepcopy(minor_idx_tmp)
        end
    end

    return minor_idx_tmp
end

function _PMinorIdx(N::Int64, sizes::Vector{Int64})
    minor_idx_dict = Dict{Int64,Vector{Tuple{Int64,Vararg{Int64}}}}()
    for k in sizes
        minor_idx_dict[k] = LOpt.get_minor_idx(N, k)
    end

    return minor_idx_dict
end

"""
    priority_central_nodes(adjacency_augment_graph::Array{<:Number}, num_nodes::Int64) 

Returns a vector of order of central nodes as 
an input to construct the graph. 
"""
function priority_central_nodes(
    adjacency_augment_graph::Array{<:Number},
    num_nodes::Int64,
)
    edge_weights_sum_list = Vector{Float64}(undef, num_nodes)
    central_nodes_list = Vector{Int64}(undef, num_nodes)

    for i in 1:num_nodes
        edge_weights_sum_list[i] = sum(adjacency_augment_graph[i, :])
    end
    central_nodes_list = sortperm(edge_weights_sum_list, rev = true)
    return central_nodes_list
end

"""
    weighted_adjacency_matrix(G::Graphs.SimpleGraphs.SimpleGraph{<:Number}, adjacency_augment_graph::Array{<:Number}, size::Int64)

Returns a weighted adjacency matrix for 
the connected part of  graph.  
"""
function weighted_adjacency_matrix(
    G::Graphs.SimpleGraphs.SimpleGraph{<:Number},
    adjacency_augment_graph::Array{<:Number},
    size::Int64,
)
    weighted_adj_matrix_size = Matrix{Float64}(undef, size, size)

    #collecting vertices connected in graph
    vertices_from_edges = Int[]
    for edge in Graphs.edges(G)
        push!(vertices_from_edges, Graphs.src(edge))
        push!(vertices_from_edges, Graphs.dst(edge))
    end
    vertices_from_edges = unique(vertices_from_edges)

    for i in 1:size
        for j in i:size
            weighted_adj_matrix_size[i, j] =
                (adjacency_augment_graph.*Graphs.adjacency_matrix(G))[
                    vertices_from_edges[i],
                    vertices_from_edges[j],
                ]
            weighted_adj_matrix_size[j, i] = weighted_adj_matrix_size[i, j]
        end
    end

    return weighted_adj_matrix_size
end

"""
    update_kopt_adjacency!(G::Graphs.SimpleGraphs.SimpleGraph{<:Number},adjacency_augment_graph::Array{<:Number},adjacency_graph_ac_tuple::Tuple{Array,Float};tol = 1E-6,)

Returns an updated adjacency_graph_ac_tuple after the kopt refinement.  
"""

function update_kopt_adjacency!(
    G::Graphs.SimpleGraphs.SimpleGraph{<:Number},
    adjacency_augment_graph::Array{<:Number},
    adjacency_graph_ac_tuple::Tuple{Array, <:Number};
    tol = 1E-6,
)
    if algebraic_connectivity(
        adjacency_augment_graph .* Matrix(Graphs.adjacency_matrix(G)),
    ) - adjacency_graph_ac_tuple[2] >= tol
        return Matrix(Graphs.adjacency_matrix(G)),
        algebraic_connectivity(
            adjacency_augment_graph .* Matrix(Graphs.adjacency_matrix(G)),
        )
    else
        return adjacency_graph_ac_tuple
    end
end

"""
    edge_combinations(n::Int64, k::Int64)

Returns all combinations possible for given `n` and `k`.
"""

function edge_combinations(
    num_edges::Int64, 
    kopt_parameter::Int64
)
    if kopt_parameter < 2 || kopt_parameter > 3
        Memento.error(_LOGGER, "kopt_parameter must be either 2 or 3")
    elseif num_edges < kopt_parameter
        Momento.error(_LOGGER, "num_edges must be greater than or equal to k")
    end

    combinations = []

    if kopt_parameter == 2
        for i in 1:num_edges-1
            for j in i+1:num_edges
                push!(combinations, (i, j))
            end
        end
    elseif kopt_parameter == 3
        for i in 1:num_edges-2
            for j in i+1:num_edges-1
                for l in j+1:num_edges
                    push!(combinations, (i, j, l))
                end
            end
        end
    end
    return combinations
end

"""
    add_multiple_edges!(G::Graphs.SimpleGraphs.SimpleGraph{<:Number},edge_set::Vector{Tuple{Int64,Int64}},)

Returns updated graph by adding edges.
"""

function add_multiple_edges!(
    G::Graphs.SimpleGraphs.SimpleGraph{<:Number},
    edge_set::Vector{Tuple{Int64,Int64}},
)
    for i in 1:length(edge_set)
        Graphs.add_edge!(G, edge_set[i][1], edge_set[i][2])
    end
    return G
end

"""
    remove_multiple_edges!(G::Graphs.SimpleGraphs.SimpleGraph{<:Number},edge_set::Vector{Tuple{Int64,Int64}},)

Returns updated graph by removing edges.
"""

function remove_multiple_edges!(
    G::Graphs.SimpleGraphs.SimpleGraph{<:Number},
    edge_set::Vector{Tuple{Int64,Int64}},
)
    for i in 1:length(edge_set)
        Graphs.rem_edge!(G, edge_set[i][1], edge_set[i][2])
    end
    return G
end
