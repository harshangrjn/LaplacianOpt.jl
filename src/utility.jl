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
        Memento.error(_LOGGER, "PSD matrix (W) cannot have complex eigenvalues")
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
function cheeger_constant(G::Matrix{<:Number}, 
                          optimizer::MOI.OptimizerWithAttributes; 
                          cheeger_ub = 1E4,
                          optimizer_log = false)

    if size(G)[1] !== size(G)[2]
        Memento.error(_LOGGER, "Input adjacency has to be a square matrix")
    end
    num_nodes = size(G)[1]
    
    result_cheeger = Dict{String, Any}("cheeger_constant" => 0)

    if !isapprox(LOpt.algebraic_connectivity(G), 0, atol = 1E-5)
        m_cheeger = JuMP.Model(optimizer)
        if !optimizer_log
            JuMP.set_silent(m_cheeger)
        end

        # Variables
        JuMP.@variable(m_cheeger, 0 <= z[1:num_nodes] <= 1, Bin)
        JuMP.@variable(m_cheeger, 0 <= Z[1:num_nodes, 1:num_nodes] <= 1)
        JuMP.@variable(m_cheeger, 0 <= cheeger <= cheeger_ub)
        JuMP.@variable(m_cheeger, cheeger_z[1:num_nodes] >= 0)

        # Constraints
        JuMP.@constraint(m_cheeger, sum(z) >= 1)
        JuMP.@constraint(m_cheeger, sum(z) <= floor(num_nodes/2))
        JuMP.@constraint(m_cheeger, sum(cheeger_z) >= sum(G[i,j] * (z[i] + z[j] -2*Z[i,j]) 
                            for i=1:(num_nodes-1), j = (i+1):num_nodes if !(isapprox(G[i,j], 0, atol=1E-5))))

        for i=1:(num_nodes-1)
            for j=1(i+1):num_nodes
                LOpt.relaxation_bilinear(m_cheeger, Z[i,j], z[i], z[j])
            end
        end
        for i=1:num_nodes
            LOpt.relaxation_bilinear(m_cheeger, cheeger_z[i], cheeger, z[i])
        end

        # Objective
        JuMP.@objective(m_cheeger, Min, cheeger)

        JuMP.optimize!(m_cheeger)

        result_cheeger["cheeger_constant"] = JuMP.objective_value(m_cheeger)
        result_cheeger["S"] = Set(findall(isapprox.(abs.(JuMP.value.(z)), 1, atol=1E-5)))
        result_cheeger["S_complement"] = Set(findall(isapprox.(abs.(JuMP.value.(z)), 0, atol=1E-5)))
        result_cheeger["solve_time"] = JuMP.solve_time(m_cheeger)
    end
    
    return result_cheeger
end