module RobinHood

using OptimizationOptimJL, ReverseDiff
using Interpolations
using NaNMath
using Base.Threads
using HDF5

include("backgrounds/background_simple.jl")
include("plots.jl")
include("actions.jl")
include("runs.jl")
include("utils.jl")

end # module