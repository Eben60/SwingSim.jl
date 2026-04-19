# """
#     Package SwingSim v$(pkgversion(SwingSim))

# Attempt of a very simplified simulation of a gymnast on a swing

# $(pathof(SwingSim))) 
# """

module SwingSim

using Modia
using GLMakie

usePlotPackage("GLMakie")


function shift_pulse(t, shift, a)
    t <= 0 && return 0.0
    t_a = sqrt(shift / a)
    if t < t_a
        0.5 * a * t^2
    elseif t < 2 * t_a
        dt = t - t_a
        shift/2 + a * dt * (t_a - dt/2)
    else
        shift
    end
end

function weightpos(l_max, shift, f, t; ϕ₀=0, pulse_dist = π/2, a=30)
    τ = 1/f
    t_red = t + 2*π*τ*ϕ₀
    t_p = mod(t_red, τ)

    # Maximum achievable shift within τ/8 total time (t_a = τ/16) is a * (τ/16)^2
    shift > a * (τ/16)^2 && error("Shift $shift cannot be achieved within time at acceleration $a")
    
    t_pulse2 = τ * (pulse_dist / (2π)) # Starts at τ/4 for default pulse_dist = π/2
    dl = shift_pulse(t_p, shift, a) - shift_pulse(t_p - t_pulse2, shift, a)
    return l_max - dl
end

function testplot(save=true; a=150)
    ts = range(0, 1, length=1000)
    ys = Float64[SwingSim.weightpos(1, 0.5, 1, t; a) for t in ts]
    fig = GLMakie.lines(ts, ys)
    GLMakie.save("plot.png", fig)
    return fig
end
export testplot

end
