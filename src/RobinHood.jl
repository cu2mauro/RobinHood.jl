module RobinHood

using StaticArrays
using OptimizationOptimJL, ReverseDiff
using Interpolations
using NaNMath
using Base.Threads
using Plots, PlotlyBase
using HDF5
using BasicBSpline

plotly()
include("backgrounds/background_simple.jl")
include("plots.jl")
include("actions.jl")
include("runs.jl")
include("utils.jl")

end # module