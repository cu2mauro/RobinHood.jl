# RobinHood

[![Build Status](https://github.com/cu2mauro/RobinHood.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/cu2mauro/RobinHood.jl/actions) [![DOI](https://zenodo.org/badge/842631903.svg)](https://zenodo.org/doi/10.5281/zenodo.13646854)

RobinHood uses numerical optimization with [Julia](https://julialang.org) to solve the Nambu-Goto holographic Wilson loops with various quivers. Specifying their rank function in a supergravity background, for each quiver it investigates confinement and screening of the theory.

## Installing

This package is not in the registry. You can install it by specifying the URL to the repository:

```julia
pkg> add https://github.com/cu2mauro/RobinHood.jl/

julia> using RobinHood
```

## Usage

To use this Julia package, first navigate to the local installation `cd("./RobinHood.jl")`. Then select the quiver in the `src/backgrounds` folder, and `include("config_file.jl")`.
