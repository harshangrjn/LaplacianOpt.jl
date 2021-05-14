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
LaplacianOpt.optimal_graph_edges() returns a vector of tuples of edges corresponding 
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
LaplacianOpt.laplacian_matrix() returns the weighted Laplacian matrix 
for an input weighted adjacency matrix of the graph. 
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
LaplacianOpt.fiedler_vector() returns the Fiedler vector or the eigenvector corresponding to the 
second smallest eigenvalue of the Laplacian matrix for an input weighted adjacency matrix of the graph. 
"""
function fiedler_vector(adjacency_matrix::Array{Float64})
    
    L_mat = LO.laplacian_matrix(adjacency_matrix)
    ac = LO.algebraic_connectivity(adjacency_matrix)

    fiedler = LA.eigvecs(L_mat)[:,2] 

    if !isapprox(fiedler' * L_mat * fiedler, ac, atol=1E-6)
        Memento.error(_LOGGER, "Evaluated eigenvector is not the Fiedler vector")
    end
    
    return fiedler

end

"""
LaplacianOpt.algebraic_connectivity() returns the algebraic connectivity or the  
second smallest eigenvalue of the Laplacian matrix, for an input weighted adjacency matrix of the graph. 
"""
function algebraic_connectivity(adjacency_matrix::Array{Float64})
    
    L_mat = LO.laplacian_matrix(adjacency_matrix)

    ac = sort(LA.eigvals(L_mat))[2] 
    
    return ac
end