@testset "relaxation_bilinear tests" begin
    LOpt.silence()

    m = JuMP.Model(glpk_optimizer)

    LB = [-1, -2.5, 0, -3, 0]
    UB = [2.5, 0, 1, 3.2, 1]
    JuMP.@variable(m, LB[i] <= x[i = 1:5] <= UB[i])
    JuMP.@variable(m, z[1:6])

    LOpt.relaxation_bilinear(m, z[1], x[1], x[2])
    LOpt.relaxation_bilinear(m, z[2], x[1], x[4])
    LOpt.relaxation_bilinear(m, z[3], x[2], x[3])
    LOpt.relaxation_bilinear(m, z[4], x[3], x[4])
    LOpt.relaxation_bilinear(m, z[5], x[3], x[5])
    LOpt.relaxation_bilinear(m, z[6], x[4], x[5])

    JuMP.@constraint(m, sum(x) >= 3.3)
    JuMP.@constraint(m, sum(x) <= 3.8)

    coefficients = [0.4852, -0.7795, 0.1788, 0.7778, 0.3418, -0.3407]
    JuMP.set_objective_function(m, sum(coefficients[i] * z[i] for i in 1:length(z)))

    obj_vals = Vector{Float64}()

    for obj_sense in [MOI.MIN_SENSE, MOI.MAX_SENSE]
        JuMP.set_objective_sense(m, obj_sense)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.OPTIMAL
        push!(obj_vals, JuMP.objective_value(m))
    end

    @test obj_vals[1] <= -9.922644 # global optimum value
    @test obj_vals[2] >= 4.908516 # global optimum value
end
