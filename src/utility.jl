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

function get_optimal_graph(num_nodes::Int, z_val::Array{Float64})

    edges = []

    for i=1:num_nodes-1
        for j=(i+1):num_nodes
           if isapprox(z_val[i,j], 1.0)
            push!(edges, (i,j))
           end
        end
    end

    return edges
end