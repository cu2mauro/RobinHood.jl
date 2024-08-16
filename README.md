# RobinHood

[![Build Status](https://github.com/cu2mauro/RobinHood.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/cu2mauro/RobinHood.jl/actions)

From various quivers of rank functions in supergravity backgrounds, use numerical optimization with Julia to solve the Nambu-Goto Wilson Loop and investigate confinement and screening


## Installing

This package is not in the registry. You can install it by specifying the URL to the repository:

```julia
pkg> add https://github.com/cu2mauro/RobinHood.jl/

julia> using RobinHood
```

## Usage

To use this Julia package, first navigate to the local installation `cd("./RobinHood.jl")`. Then select the quiver in the `src/backgrounds` folder, and `include("config_file.jl")`.
