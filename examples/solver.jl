"""
    get_solver(params)
    This function returns the JuMP optimizer with its appropriate attributes, based on 
        user-defined inputs in `params` dictionary.  
"""
function get_solver(params::Dict{String,Any})
    optimizers_list = ["cplex", "cbc", "ipopt", "juniper", "alpine", "glpk", "gurobi"]

    # Optimizer
    if !("optimizer" in keys(params))
        error("Input a valid MIP optimizer")
    end

    if !(params["optimizer"] in optimizers_list)
        error("Specified optimizer does not belong in the pre-defined list. Add your optimizer separately with it's attributes")
    end

    if "presolve" in keys(params)
        presolve = params["presolve"]
    else
        # default value
        presolve = true
    end

    if "optimizer_log" in keys(params)
        optimizer_log = params["optimizer_log"]
    else
        # default value
        optimizer_log = true
    end

    # Mixed-integer programming optimizers
    if params["optimizer"] == "cplex"    # commercial 
       return get_cplex(presolve, optimizer_log) 
       
    elseif params["optimizer"] == "gurobi"    # commercial 
        return get_gurobi(presolve, optimizer_log) 

    elseif params["optimizer"] == "cbc"  # open-source
       return get_cbc(presolve, optimizer_log)

    elseif params["optimizer"] == "glpk"  # open-source
        return get_glpk(presolve, optimizer_log)
    
    # Local mixed-integer nonlinear programming optimizers
    elseif params["optimizer"] == "ipopt"    # open-source 
       return get_ipopt(presolve, optimizer_log)
       
    elseif params["optimizer"] == "juniper"  # open-source 
       return get_juniper(presolve, optimizer_log)
    
    # Global NLP/MINLP optimizer
    elseif params["optimizer"] == "alpine"   # open-source 
         alpine = JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            "nlp_solver" => get_ipopt(presolve, optimizer_log),
                                            "minlp_solver" => get_juniper(presolve, optimizer_log),  
                                            "mip_solver" => get_cplex(presolve, optimizer_log),
                                            "presolve_bt" => false,
                                            "presolve_max_iter" => 10,
                                            "presolve_bp" => false,
                                            "disc_ratio" => 10)
        return alpine
    end
end

function get_cplex(presolve::Bool, optimizer_log::Bool)
     cplex = JuMP.optimizer_with_attributes(CPLEX.Optimizer, 
                                       MOI.Silent() => !optimizer_log, 
                                       "CPX_PARAM_PREIND" => presolve) 
    return cplex 
end

function get_gurobi(presolve::Bool, optimizer_log::Bool)
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer, 
                                      MOI.Silent() => !optimizer_log, 
                                      "Presolve" => presolve) 
   return gurobi
end

function get_ipopt(presolve::Bool, optimizer_log::Bool)
     ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
                                       MOI.Silent() => !optimizer_log, 
                                       "sb" => "yes", 
                                       "max_iter" => 1E4)
    return ipopt 
end

function get_cbc(presolve::Bool, optimizer_log::Bool)
     cbc = JuMP.optimizer_with_attributes(Cbc.Optimizer, 
                                     MOI.Silent() => !optimizer_log)
    return cbc 
end

function get_juniper(presolve::Bool, optimizer_log::Bool)
     juniper = JuMP.optimizer_with_attributes(Juniper.Optimizer, 
                                         MOI.Silent() => !optimizer_log, 
                                         "mip_solver" => get_cplex(presolve, optimizer_log), 
                                         "nl_solver" => get_ipopt(presolve, optimizer_log))
    return juniper
end

function get_glpk(presolve::Bool, optimizer_log::Bool)
    glpk = JuMP.optimizer_with_attributes(GLPK.Optimizer, 
                                     MOI.Silent() => !optimizer_log)
   return glpk
end
