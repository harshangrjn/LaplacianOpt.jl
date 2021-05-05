module LaplacianOpt

import JuMP
import LinearAlgebra
import Random
import Memento
import MathOptInterface
import JSON
import LightGraphs
import TikzGraphs
import TikzPictures

const MOI = MathOptInterface
const LA = LinearAlgebra
const LO = LaplacianOpt
const LG = LightGraphs

# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(LaplacianOpt)`
# NOTE: If this line is not included then the precompiled `LaplacianOpt._LOGGER` won't be registered at runtime.
__init__() = Memento.register(_LOGGER)

"Suppresses information and warning messages output by LaplacianOpt, for fine grained control use the Memento package"
function silence()
    Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(LaplacianOpt), "error")
end

"allows the user to set the logging level without the need to add Memento"
function logger_config!(level)
    Memento.config!(Memento.getlogger("LaplacianOpt"), level)
end

include("data.jl")
include("types.jl")
include("utility.jl")
include("relaxations.jl")
include("lopt_model.jl")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("solution.jl")
include("log.jl")

end # module
