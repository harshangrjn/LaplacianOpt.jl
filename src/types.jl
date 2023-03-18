export LaplacianOptModel

"""
    LaplacianOptModelOptions
The composite mutable struct, `LaplacianOptModelOptions`, holds various optimization model options for enhancements 
with defualt options set to the values provided by `get_default_options` function.
"""
mutable struct LaplacianOptModelOptions
    solution_type::String
    formulation_type::String

    eigen_cuts_full::Bool
    eigen_cuts_2minors::Bool
    eigen_cuts_3minors::Bool
    projected_eigen_cuts::Bool
    soc_linearized_cuts::Bool
    topology_multi_commodity::Bool
    topology_flow_cuts::Bool
    cheeger_cuts::Bool
    cheeger_cuts_factor::Float64

    sdp_relaxation::Bool

    kopt_parameter::Int
    best_lower_bound::Float64
    best_incumbent::Union{Vector{Tuple{Int64,Int64}},Nothing}

    tol_zero::Float64
    tol_psd::Float64
    time_limit::Float64
    relax_integrality::Bool
    lazycuts_logging::Bool
end

"""
    get_default_options()
This function returns the default options for building the struct `LaplacianOptModelOptions`.
"""
function get_default_options()
    solution_type = "optimal" # optimal, heuristic
    formulation_type = "max_λ2"  # max_λ2, max_span_tree

    eigen_cuts_full = true   # true, false     
    eigen_cuts_2minors = true   # true, false
    eigen_cuts_3minors = false  # true, false
    projected_eigen_cuts = true  # true, false
    soc_linearized_cuts = false  # true, false
    topology_multi_commodity = false  # true, false
    topology_flow_cuts = true   # true, false
    cheeger_cuts = false  # true, false
    cheeger_cuts_factor = 0.5    # > 0.5 may make cheeger_cuts a heuristic cut

    sdp_relaxation = false  # true, false

    kopt_parameter = 2      # Integer value >= 1
    best_lower_bound = 0      # Best known feasible solution's objective
    best_incumbent = nothing

    tol_zero = 1E-6   # > 0 value
    tol_psd = 1E-6   # > 0 value (tolerance to verify PSD-ness of a matrix)
    time_limit = 10800  # > 0 value (seconds)
    relax_integrality = false  # true, false
    lazycuts_logging = false  # true, false

    return LaplacianOptModelOptions(
        solution_type,
        formulation_type,
        eigen_cuts_full,
        eigen_cuts_2minors,
        eigen_cuts_3minors,
        projected_eigen_cuts,
        soc_linearized_cuts,
        topology_multi_commodity,
        topology_flow_cuts,
        cheeger_cuts,
        cheeger_cuts_factor,
        sdp_relaxation,
        kopt_parameter,
        best_lower_bound,
        best_incumbent,
        tol_zero,
        tol_psd,
        time_limit,
        relax_integrality,
        lazycuts_logging,
    )
end

"""
    LaplacianOptModel

The composite mutable struct, `LaplacianOptModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct LaplacianOptModel
    data::Dict{String,Any}
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    options::LaplacianOptModelOptions
    result::Dict{String,Any}

    "Contructor for struct `LaplacianOptModel`"
    function LaplacianOptModel(data::Dict{String,Any})
        data = data
        model = JuMP.Model()
        variables = Dict{Symbol,Any}()
        options = LOpt.get_default_options()
        result = Dict{String,Any}("cheeger_cuts_separation_time" => 0)
        lom = new(data, model, variables, options, result)

        return lom
    end
end

"""
    GraphData

The composite mutable struct, `GraphData`, holds matrices of adjacency matrix, laplacian matrix,
fiedler vector and algebraic connectivity.
"""
mutable struct GraphData
    adjacency::Array{<:Number}
    laplacian::Array{<:Number}
    fiedler::Vector{<:Number}
    ac::Number

    function GraphData(adjacency::Array{<:Number})
        adj_matrix = adjacency
        laplacian = LOpt.laplacian_matrix(adjacency)
        fiedler = LOpt.fiedler_vector(adjacency)
        ac = LOpt.algebraic_connectivity(adjacency)

        graph_data = new(adj_matrix, laplacian, fiedler, ac)

        return graph_data
    end
end
