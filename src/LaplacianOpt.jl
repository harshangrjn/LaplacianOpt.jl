module LaplacianOpt

import JuMP
import LinearAlgebra
import Logging
import MathOptInterface
import JSON
import Graphs
import TikzGraphs
import TikzPictures
import Pkg

const MOI = MathOptInterface
const LA = LinearAlgebra
const LOpt = LaplacianOpt

# Create a package-local logger so LaplacianOpt logging can be configured
# without changing the user's global logger.
const _LOGGER = Ref{Logging.ConsoleLogger}()

function __init__()
    logger_config!("info")
    return
end

"Suppresses information and warning messages output by LaplacianOpt."
function silence()
    _log_if_level(
        () -> @info(
            "Suppressing information and warning messages for the rest of this session."
        ),
        Logging.Info,
    )
    return logger_config!("error")
end

"allows the user to set the logging level"
function logger_config!(level::Logging.LogLevel)
    _LOGGER[] = Logging.ConsoleLogger(stdout, level; meta_formatter = _meta_formatter)
    return
end

function logger_config!(level::String)
    return getfield(Logging, level |> titlecase |> Symbol) |> logger_config!
end

function _meta_formatter(level::Logging.LogLevel, _module, args...)
    return Logging.default_logcolor(level), "$(_module) | $level]:", ""
end

function _log_if_level(f, level, logger = _LOGGER[])
    if level >= Logging.min_enabled_level(logger)
        Logging.with_logger(f, logger)
    end
    return
end

macro _error(msg)
    return quote
        $_log_if_level(() -> @error($msg), $(Logging.Error))
        error($msg)
    end |> esc
end

macro _warn(msg)
    return :($_log_if_level(() -> @warn($msg), $(Logging.Warn))) |> esc
end

macro _debug(msg)
    return :($_log_if_level(() -> @debug($msg), $(Logging.Debug))) |> esc
end

macro _info(msg)
    return :($_log_if_level(() -> @info($msg), $(Logging.Info))) |> esc
end

include("data.jl")
include("types.jl")
include("utility.jl")
include("relaxations.jl")
include("lopt_model.jl")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("heuristics.jl")
include("solution.jl")
include("log.jl")

end # module
