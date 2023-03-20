#=====================================#
# MIP solvers (commercial, but fast)  #
#=====================================#
function get_gurobi()
    GRB_ENV = Gurobi.Env()
    return optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV),
        MOI.Silent() => true,
        # "MIPFocus" => 3, # Focus on optimality over feasibility 
        "Presolve" => 1,
    )
end

function get_cplex()
    return JuMP.optimizer_with_attributes(
        CPLEX.Optimizer,
        MOI.Silent() => true,
        # "CPX_PARAM_EPGAP" => 1E-6,
        # "CPX_PARAM_EPOPT" => 1E-8,
        # "CPX_PARAM_MIPEMPHASIS" => 1, # Focus on optimality over feasibility 
        "CPX_PARAM_PREIND" => 1,
    )
end

#======================================#
# MIP solvers (open-source, but slow)  #
#======================================#

function get_glpk()
    return JuMP.optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => false)
end

# https://github.com/jump-dev/HiGHS.jl
function get_highs()
    return JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "presolve" => "on",
        "log_to_console" => true,
    )
end

#========================================================#
# Continuous nonlinear programming solver (open-source)  #
#========================================================#

function get_ipopt()
    return JuMP.optimizer_with_attributes(
        Ipopt.Optimizer,
        MOI.Silent() => true,
        "sb" => "yes",
        "max_iter" => Int(1E4),
    )
end

#=================================================================#
# Local mixed-integer nonlinear programming solver (open-source)  #
#=================================================================#
function get_juniper()
    return JuMP.optimizer_with_attributes(
        Juniper.Optimizer,
        MOI.Silent() => false,
        "mip_solver" => get_gurobi(),
        "nl_solver" => get_ipopt(),
    )
end
