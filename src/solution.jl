function build_LOModel_result(lom::LaplacianOptModel, solve_time::Number)
    # try-catch is needed until solvers reliably support ResultCount()
    result_count = 1
    try
        result_count = JuMP.result_count(lom.model)
    catch
        Memento.warn(
            _LOGGER,
            "the given optimizer does not provide the ResultCount() attribute, assuming the solver returned a solution which may be incorrect.",
        )
    end

    solution = Dict{String,Any}()

    if result_count > 0
        solution = LOpt.build_LOModel_solution(lom)
    else
        Memento.warn(
            _LOGGER,
            "LaplacianOpt model has no results - solution cannot be built",
        )
    end

    result = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(lom.model),
        "termination_status" => JuMP.termination_status(lom.model),
        "primal_status" => JuMP.primal_status(lom.model),
        "objective" => LOpt.get_objective_value(lom.model),
        "objective_ub" => LOpt.get_objective_bound(lom.model),
        "solve_time" => solve_time,
        "solution" => solution,
        "solution_type" => lom.options.solution_type,
        "adjacency_base_graph" => lom.data["adjacency_base_graph"],
        "adjacency_augment_graph" => lom.data["adjacency_augment_graph"],
    )

    status = LOpt.optimality_certificate_MISDP(lom, result)
    if status in [true, false]
        result["optimality_certificate_MISDP"] = status
    end

    return result
end

""
function get_objective_value(model::JuMP.Model)
    obj_val = NaN

    try
        obj_val = JuMP.objective_value(model)
    catch
        Memento.warn(
            _LOGGER,
            "Objective value is unbounded. Problem may be infeasible or not constrained properly",
        )
    end

    return obj_val
end

""
function get_objective_bound(model::JuMP.Model)
    obj_lb = -Inf

    try
        obj_lb = JuMP.objective_bound(model)
    catch
    end

    return obj_lb
end

function build_LOModel_solution(lom::LaplacianOptModel)
    solution = Dict{String,Any}()

    for i in keys(lom.variables)
        solution[String(i)] = JuMP.value.(lom.variables[i])
    end

    return solution
end

"""
    optimality_certificate_MISDP(lom::LaplacianOptModel, result::Dict{String,Any})

Given the LaplacianOpt model (`lom`) and the results dictionary, this function returns 
a boolean if the integer solution satisfies the optimality certificate for the MISDP problem. 
"""
function optimality_certificate_MISDP(lom::LaplacianOptModel, result::Dict{String,Any})
    if ("z_var" in keys(result["solution"])) &&
       !(lom.options.formulation_type == "max_span_tree")
        adjacency_full_graph = lom.data["adjacency_augment_graph"]
        (lom.data["num_edges_existing"] > 0) &&
            (adjacency_full_graph += lom.data["adjacency_base_graph"])
        ac = LOpt.algebraic_connectivity(
            abs.(result["solution"]["z_var"]) .* adjacency_full_graph,
        )
        if isapprox(ac, result["objective"], atol = 1E-5)
            return true
        else
            return false
        end
    end
end

function build_LOModel_heuristic_result(lom::LaplacianOptModel, solve_time::Number)

    solution = Dict{String,Any}()

    result = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(lom.model),
        # "objective" => LOpt.get_objective_value(lom.model),
        "solve_time" => solve_time,
        "solution" => solution,
        "solution_type" => lom.options.solution_type,
        "adjacency_base_graph" => lom.data["adjacency_base_graph"],
        "adjacency_augment_graph" => lom.data["adjacency_augment_graph"],
    )

    return result
end