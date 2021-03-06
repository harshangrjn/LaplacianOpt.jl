export LaplacianOptModel

"""
    LaplacianOptModel

The composite mutable struct, `LaplacianOptModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct LaplacianOptModel 
    data      :: Dict{String,Any}
    model     :: JuMP.Model
    variables :: Dict{Symbol,Any}
    result    :: Dict{String,Any}

    "Contructor for struct `LaplacianOptModel`"
    function LaplacianOptModel(data::Dict{String,Any})
        
        data = data
        model = JuMP.Model() 
        variables = Dict{Symbol,Any}()
        result = Dict{String,Any}()
        lom = new(data, model, variables, result)

        return lom
    end

end

"""
    GraphData

The composite mutable struct, `GraphData`, holds matrices of adjacency matrix, laplacian matrix,
fiedler vector and algebraic connectivity.
"""
mutable struct GraphData
    adjacency ::Array{<:Number}
    laplacian ::Array{<:Number}
    fiedler   ::Vector{<:Number}
    ac        ::Number
    
    function GraphData(adjacency::Array{<:Number})
        adj_matrix = adjacency
        laplacian = LOpt.laplacian_matrix(adjacency)
        fiedler = LOpt.fiedler_vector(adjacency)
        ac = LOpt.algebraic_connectivity(adjacency)

        graph_data = new(adj_matrix, laplacian, fiedler, ac)

        return graph_data
        
    end
end
