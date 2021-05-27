export LaplacianOptModel

"""
The composite mutable struct, `LaplacianOptModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct LaplacianOptModel 
    data::Dict{String,Any}
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    result::Dict{String,Any}

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

