"""
    round_zeros_ones!(v::Array{<:Number}; tol = 1E-6)

Given a vector of numbers, this function updates the input vector by rounding the values closest to 0, 1 and -1. 
"""
function round_zeros_ones!(v::Array{<:Number}; tol = 1E-6)
    v[abs.(v) .<= tol] .= 0
    near_ones = findall(((1 - tol) .<= abs.(v) .<= (1 + tol)))
    return v[near_ones] .= sign.(v[near_ones])
end

"""
    optimal_graph_edges(adjacency_matrix::Array{<:Number}) 

Returns a vector of tuples of edges corresponding 
to an input adjacency matrix of the graph. 
"""
function optimal_graph_edges(adjacency_matrix::Array{<:Number})
    num_nodes = size(adjacency_matrix, 1)

    !LA.issymmetric(adjacency_matrix) &&
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")

    any(i -> !isapprox(adjacency_matrix[i, i], 0), 1:num_nodes) &&
        Memento.error(_LOGGER, "Input adjacency matrix cannot have self loops")

    edges = Tuple{Int64,Int64}[]
    for i in 1:(num_nodes-1), j in (i+1):num_nodes
        !isapprox(adjacency_matrix[i, j], 0) && push!(edges, (i, j))
    end

    return edges
end

"""
    laplacian_matrix(adjacency_matrix::Array{<:Number}) 

Given a weighted adjacency matrix as an input, this function returns the 
weighted Laplacian matrix of the graph. 
"""
function laplacian_matrix(adjacency_matrix::Array{<:Number})
    !LA.issymmetric(adjacency_matrix) &&
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")

    num_nodes = size(adjacency_matrix, 1)
    laplacian_matrix = zeros(Float64, num_nodes, num_nodes)

    # Check for self loops
    for i in 1:num_nodes
        !isapprox(adjacency_matrix[i, i], 0) &&
            Memento.error(_LOGGER, "Input adjacency matrix cannot have self loops")
    end

    # Negative weights
    any(x -> x < -1E-6, adjacency_matrix) &&
        Memento.error(_LOGGER, "Input adjacency matrix cannot have negative weights")

    # Degree matrix
    for i in 1:num_nodes
        laplacian_matrix[i, i] = sum(adjacency_matrix[i, :])
    end

    # Off-diagonal elements
    for i in 1:num_nodes, j in (i+1):num_nodes
        if !isapprox(adjacency_matrix[i, j], 0)
            laplacian_matrix[i, j] = -adjacency_matrix[i, j]
            laplacian_matrix[j, i] = -adjacency_matrix[i, j]
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
    count = 0
    for i in cutset_f, j in cutset_t
        if !isapprox(adjacency[i, j], 0, atol = 1E-6)
            count += 1
            count >= num_edges && return true
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

    # Check for complex eigenvalues
    if eltype(W_eigvals) <: Complex &&
       !isapprox(imag.(W_eigvals), zeros(length(W_eigvals)), atol = 1E-6)
        Memento.error(_LOGGER, "PSD matrix (W) cannot have complex eigenvalues")
    end

    # If minimum eigenvalue is negative (beyond tolerance)
    min_eigval = LA.eigmin(W)
    if min_eigval <= -tol
        # Check if the first eigenvalue is the minimum
        if isapprox(W_eigvals[1], min_eigval, atol = 1E-6)
            violated_eigen_vec = LA.eigvecs(W)[:, 1]
            LOpt.round_zeros_ones!(violated_eigen_vec)
            return violated_eigen_vec
        else
            Memento.warn(
                _LOGGER,
                "Eigen cut corresponding to the negative eigenvalue could not be evaluated",
            )
        end
    end
    return nothing
end

"""
    cheeger_constant(G::Matrix{<:Number}, 
                     optimizer::MOI.OptimizerWithAttributes; 
                     cheeger_ub = 1E6,
                     optimizer_log = false)
    
Given a symmetric, weighted adjacency matrix `G` of a connected graph, and a mixed-integer linear optimizer, 
this function returns the Cheeger's constant (or isoperimetric number of the graph) and the associated two partitions 
(`S` and `S_complement`) of the graph. 

The exact version of the formulation implemented here is published in:
Somisetty, N., Nagarajan, H. and Darbha, S., 2024. "Spectral Graph Theoretic Methods for Enhancing Network Robustness 
in Robot Localization," IEEE 63rd Conference on Decision and Control (CDC), 2024. 
Link: https://doi.org/10.1109/CDC56724.2024.10886863

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
    get_minor_idx(vector::Vector{Int64}, size::Int64)

Given a vector and the desired size of tuples, this 
recursive function outputs all possible tuples of the 
specified size from the elements of the vector.  
"""
function get_minor_idx(vector::Vector{Int64}, size::Int64)
    if size == 1
        return [(x,) for x in vector]
    elseif size == length(vector)
        return [tuple(vector...)]
    else
        tuples = []
        for i in 1:length(vector)
            element = vector[i]
            sub_tuples = LOpt.get_minor_idx(vector[(i+1):end], size - 1)
            for sub_tuple in sub_tuples
                push!(tuples, (element, sub_tuple...))
            end
        end

        return tuples
    end
end

function _PMinorIdx(
    N::Int64,
    sizes::Vector{Int64},
    minors_on_augment_edges::Bool,
    data::Dict{String,Any},
)
    minor_idx_dict = Dict{Int64,Vector{Tuple{Int64,Vararg{Int64}}}}()

    isempty(sizes) && return minor_idx_dict

    max_size = maximum(sizes)
    max_size > N &&
        Memento.warn(_LOGGER, "Detected maximum eigen-cut size ($(max_size)) > num_nodes")

    for k in unique(sizes)
        (k < 2 || k > N) && continue

        vertices = collect(1:N)
        if minors_on_augment_edges && k < N
            vertices = LOpt._vertices_with_augment_edges(data)
        end

        minor_idx_dict[k] = LOpt.get_minor_idx(vertices, k)
    end

    return minor_idx_dict
end

function _vertices_with_augment_edges(data::Dict{String,Any})
    v = Set{Int64}()
    adj = data["adjacency_augment_graph"]
    for i in 1:(data["num_nodes"]-1), j in (i+1):data["num_nodes"]
        if !isapprox(adj[i, j], 0, atol = 1E-6)
            push!(v, i, j)
        end
    end
    return collect(v)
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
    edge_weights = [sum(adjacency_augment_graph[i, :]) for i in 1:num_nodes]
    return sortperm(edge_weights, rev = true)
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
    # Get unique vertices from edges
    vertices_from_edges =
        unique(vcat([[Graphs.src(e), Graphs.dst(e)] for e in Graphs.edges(G)]...))

    # Pre-compute the product matrix
    product_matrix = adjacency_augment_graph .* Graphs.adjacency_matrix(G)

    # Initialize result matrix
    weighted_adj_matrix_size = zeros(Float64, size, size)

    # Fill upper triangular part and mirror to lower triangular
    for i in 1:size, j in i:size
        val = product_matrix[vertices_from_edges[i], vertices_from_edges[j]]
        weighted_adj_matrix_size[i, j] = val
        weighted_adj_matrix_size[j, i] = val
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
    adjacency_graph_ac_tuple::Tuple{Array,<:Number};
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

function edge_combinations(num_edges::Int64, kopt_parameter::Int64)
    if kopt_parameter < 2 || kopt_parameter > 3
        Memento.error(_LOGGER, "kopt_parameter must be either 2 or 3")
    elseif num_edges < kopt_parameter
        Memento.error(_LOGGER, "num_edges must be greater than or equal to k")
    end

    if kopt_parameter == 2
        return [(i, j) for i in 1:(num_edges-1) for j in (i+1):num_edges]
    else # kopt_parameter == 3
        return [
            (i, j, l) for i in 1:(num_edges-2) for j in (i+1):(num_edges-1) for
            l in (j+1):num_edges
        ]
    end
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
