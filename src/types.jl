export LaplacianOptModel

"""
    LaplacianOptModelOptions
The composite mutable struct, `LaplacianOptModelOptions`, holds various optimization 
model options for enhancements with defualt options set to the values provided by 
`get_default_options` function.
"""
mutable struct LaplacianOptModelOptions
    solution_type::String
    formulation_type::String

    eigen_cuts_sizes::Vector{Int64}
    projected_eigen_cuts::Bool
    soc_linearized_cuts::Bool
    minors_on_augment_edges::Bool
    topology_multi_commodity::Bool
    topology_flow_cuts::Bool
    max_cycle_length::Int64
    cheeger_cuts::Bool
    cheeger_cuts_factor::Float64

    sdp_relaxation::Bool

    kopt_parameter::Int64
    num_central_nodes_kopt::Int64
    num_swaps_bound_kopt::Int64
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

    eigen_cuts_sizes = [] #vector of integer values from 1 up to num_nodes
    projected_eigen_cuts = true  # true, false
    soc_linearized_cuts = false  # true, false
    minors_on_augment_edges = true # true, false (applies to only principal minors of sizes < num_nodes)
    topology_multi_commodity = false  # true, false
    topology_flow_cuts = true   # true, false
    max_cycle_length = 3 # >= 0 and <= (num_nodes-1) integer value 
    cheeger_cuts = false  # true, false
    cheeger_cuts_factor = 0.5    # > 0.5 may make cheeger_cuts a heuristic cut

    sdp_relaxation = false  # true, false

    kopt_parameter = 2      # 1 <= Integer value <= 3
    num_central_nodes_kopt = 5 # 1 <= Integer value <= size of instance
    num_swaps_bound_kopt = 1E4  # Integer value  
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
        eigen_cuts_sizes,
        projected_eigen_cuts,
        soc_linearized_cuts,
        minors_on_augment_edges,
        topology_multi_commodity,
        topology_flow_cuts,
        max_cycle_length,
        cheeger_cuts,
        cheeger_cuts_factor,
        sdp_relaxation,
        kopt_parameter,
        num_central_nodes_kopt,
        num_swaps_bound_kopt,
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
    minor_idx_dict::Dict{Int64,Vector{Tuple{Int64,Vararg{Int64}}}}
    result::Dict{String,Any}

    "Contructor for struct `LaplacianOptModel`"
    function LaplacianOptModel(data::Dict{String,Any})
        data = data
        model = JuMP.Model()
        variables = Dict{Symbol,Any}()
        options = LOpt.get_default_options()
        minor_idx_dict = Dict{Int64,Vector{Tuple{Int64,Vararg{Int64}}}}()
        result = Dict{String,Any}("cheeger_cuts_separation_time" => 0)
        lom = new(data, model, variables, options, minor_idx_dict, result)

        return lom
    end
end

"""
    GraphData

The composite mutable struct, `GraphData`, holds matrices of adjacency matrix, 
laplacian matrix, fiedler vector and algebraic connectivity.
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
