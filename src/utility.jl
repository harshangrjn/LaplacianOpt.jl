function get_rounded_zeros_and_ones!(v::Array{Float64}, tol_zero::Float64)
   
    for i in findall(abs.(v) .<= tol_zero)
        v[i] = 0
    end

    for i in findall(((1-tol_zero) .<= abs.(v) .<= (1+tol_zero)))
        if v[i] > 0 
            v[i] = 1
        elseif v[i] < 0 
            v[i] = -1
        end
    end
end 

"""
    optimal_graph_edges(adjacency_matrix::Array{Float64}) 

Returns a vector of tuples of edges corresponding 
to an input adjacency matrix of the graph. 
"""
function optimal_graph_edges(adjacency_matrix::Array{Float64})

    edges = Vector{Tuple{Int64, Int64}}()
    
    num_nodes = size(adjacency_matrix)[1]

    if !LA.issymmetric(adjacency_matrix)
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")
    end

    for i=1:num_nodes
        if !isapprox(adjacency_matrix[i,i], 0)
            Memento.error(_LOGGER, "Input adjacency matrix cannot have self loops")
        end
    end

    for i=1:(num_nodes-1)
        for j=(i+1):num_nodes
            if !isapprox(adjacency_matrix[i,j], 0)
                push!(edges, (i,j))
            end
        end
    end

    return edges
end

"""
    laplacian_matrix(adjacency_matrix::Array{Float64}) 

Returns the weighted Laplacian matrix for an input weighted adjacency matrix of the graph. 
"""
function laplacian_matrix(adjacency_matrix::Array{Float64})

    if !LA.issymmetric(adjacency_matrix)
        Memento.error(_LOGGER, "Input adjacency matrix is asymmetric")
    end

    num_nodes = size(adjacency_matrix)[1]

    laplacian_matrix = zeros(Float64, num_nodes, num_nodes)

    for i = 1:num_nodes
        for j = i:num_nodes
            
            if i == j    
                if !isapprox(adjacency_matrix[i,j], 0)
                    Memento.error(_LOGGER, "Input adjacency matrix cannot have self loops")
                end

                laplacian_matrix[i,j] = sum(adjacency_matrix[i,:])
            else 
                if adjacency_matrix[i,j] <= -1E-6
                    Memento.error(_LOGGER, "Input adjacency matrix cannot have negative weights")
                end

                if !isapprox(adjacency_matrix[i,j], 0)
                    laplacian_matrix[i,j] = -adjacency_matrix[i,j]
                    laplacian_matrix[j,i] = -adjacency_matrix[i,j]
                end
            end

        end
    end

    return laplacian_matrix
end

"""
    fiedler_vector(adjacency_matrix::Array{Float64})

Returns the Fiedler vector or the eigenvector corresponding to the 
second smallest eigenvalue of the Laplacian matrix for an input weighted adjacency matrix of the graph. 
"""
function fiedler_vector(adjacency_matrix::Array{Float64})
    
    L_mat = LOpt.laplacian_matrix(adjacency_matrix)
    ac    = LOpt.algebraic_connectivity(adjacency_matrix)

    fiedler = LA.eigvecs(L_mat)[:,2] 

    if !isapprox(fiedler' * L_mat * fiedler, ac, atol=1E-6)
        Memento.error(_LOGGER, "Evaluated eigenvector is not the Fiedler vector")
    end
    
    return fiedler

end

"""
    algebraic_connectivity(adjacency_matrix::Array{Float64}) 
    
Returns the algebraic connectivity or the  
second smallest eigenvalue of the Laplacian matrix, for an input weighted adjacency matrix of the graph. 
"""
function algebraic_connectivity(adjacency_matrix::Array{Float64})
    
    L_mat = LOpt.laplacian_matrix(adjacency_matrix)

    ac = sort(LA.eigvals(L_mat))[2] 
    
    return ac
end

function _is_flow_cut_valid(cutset_f::Vector{Int64}, cutset_t::Vector{Int64}, adjacency::Array{Float64})
    for i in cutset_f
        for j in cutset_t
            if !(isapprox(adjacency[i,j], 0, atol = 1E-6))
                return true
            end
        end
    end

    return false
end

"""
    _violated_eigen_vector(W::Array{Float64}; tol = 1E-6)
    
Given a square symmetric matrix, this function returns an eigen vector corresponding to it's most 
violated eigen value w.r.t positive semi-definiteness of the input matrix.  
"""
function _violated_eigen_vector(W::Array{Float64}; tol = 1E-6)
    W_eigvals = LA.eigvals(W)

    if typeof(W_eigvals) != Vector{Float64}
        Memento.error(_LOGGER, "PSD matrix (W) cannot have complex eigenvalues")
    end

    if LA.eigmin(W) <= -tol
        
        if W_eigvals[1] == LA.eigmin(W)
            violated_eigen_vec = LA.eigvecs(W)[:,1]
            LOpt.get_rounded_zeros_and_ones!(violated_eigen_vec, 1E-6)
            return violated_eigen_vec
        else
            Memento.warn(_LOGGER, "Eigen cut corresponding to the negative eigenvalue could not be evaluated")
        end
    else
        return
    end
end