module RobinHood

using OptimizationOptimJL, ReverseDiff
using Interpolations
using NaNMath
using Base.Threads
using Plots
using HDF5
using BasicBSpline, StaticArrays

include("backgrounds/background_simple.jl")
include("plots.jl")
include("actions.jl")
include("runs.jl")
include("utils.jl")

end # module