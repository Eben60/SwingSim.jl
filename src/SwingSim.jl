"""
    Package SwingSim v$(pkgversion(SwingSim))

Attempt of a very simplified simulation of a gymnast on a swing

$(pathof(SwingSim))) 
"""
module SwingSim
using Plots, OrdinaryDiffEq

const g = 9.81

include("weightposition.jl")

end
