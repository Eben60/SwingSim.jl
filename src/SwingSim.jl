"""
    Package SwingSim v$(pkgversion(SwingSim))

Attempt of a very simplified simulation of a gymnast on a swing

$(pathof(SwingSim))) 
"""
module SwingSim

using OrdinaryDiffEq
using ModelingToolkit
using Plots

const g = 9.81

include("weightposition.jl")
export weightspos_testplot

include("Pendulum_ODE.jl")
export solve_pendulum, plot_pendulum_polar_ODE

include("Pendulum_Model.jl")
export pendulum_xy_model, plot_pendulum_xy_model

end
