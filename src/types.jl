export LaplacianOptModel

"""
The composite mutable struct, `LaplacianOptModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct LaplacianOptModel 
    data::Dict{String,Any}
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    #constraints::Dict{Symbol,Any}
    result::Dict{String,Any}
end

"Contructor for struct `LaplacianOptModel`"
function LaplacianOptModel(data::Dict{Symbol,Any};)
    lom = new(data, JuMP.Model(), Dict{Symbol,Any}(), Dict{Symbol,Any}())
    lom.data = data
    lom.model = JuMP.Model() 
    lom.variables = Dict{Symbol,Any}()
    lom.result = Dict{String,Any}()
    return lom
end

