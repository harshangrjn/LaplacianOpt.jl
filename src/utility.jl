function get_rounded_zeros_and_ones!(v::Vector{Float64}, tol_zero::Float64)
    
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