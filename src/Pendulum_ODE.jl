function pendulum_ode(dθ, θ, p, t)
    l_max, shift, τ, ϕ₀, pulse_dist, g_no, Q, L0, span = p
    
    l = weightpos(t, l_max, shift, τ, span; ϕ₀, pulse_dist, g_no)
    dl = weightpos_prim(t, l_max, shift, τ, span; ϕ₀, pulse_dist, g_no)
    
    damping = (2 * dl / l) + (1 / Q) * sqrt(g / L0)
    
    return -damping * dθ - (g / l) * sin(θ)
end

"""
    solve_pendulum(θ₀, ω₀, tspan, l_max, shift, τ; kwargs...)

Constructs and solves the SecondOrderODEProblem for the damped pendulum with variable length.
"""
function solve_pendulum(θ₀, ω₀, tspan, l_max, shift, τ; 
                        ϕ₀=0, pulse_dist=π/2, g_no=2, 
                        Q=Inf, L0=l_max, 
                        span=nothing,
                        solver=Tsit5(), kwargs...)
    
    p = (l_max, shift, τ, ϕ₀, pulse_dist, g_no, Q, L0, span)
    prob = SecondOrderODEProblem(pendulum_ode, ω₀, θ₀, tspan, p)
    
    return solve(prob, solver; kwargs...)
end

"""
    solve_testplot(; kwargs...)

Tests solve_pendulum with default parameters and plots both angle θ and length L over time.
"""
function solve_testplot(; θ₀=-π/4, ω₀=0.0, l_max=3, τ=3.5, shift=0.2, g_no=1, ϕ₀=0, pulse_dist=π/2, Q=10, tplot=(0.0, 3*τ), solver=BS3(), kwargs...)

    p_start, p_end = tplot
    tspan = (0, p_end)

    sol = solve_pendulum(θ₀, ω₀, tspan, l_max, shift, τ; 
                         ϕ₀, pulse_dist, g_no, Q, L0=l_max, solver, kwargs...)
    
    ts = range(p_start, p_end, length=1000)
    # sol(t)[2] is the position (angle θ), sol(t)[1] is the velocity (ω)
    θs = [sol(t)[2] for t in ts]
    ls = [weightpos(t, l_max, shift, τ; ϕ₀, pulse_dist, g_no) for t in ts]
    
    theta0_deg = round(Int, rad2deg(θ₀))
    phi0_deg = round(Int, rad2deg(ϕ₀))

    p1 = plot(ts, θs, label="θ (angle)", xlabel="t", ylabel="θ (rad)", lw=2)
    annotate!(p1, (ts[1] + ts[end])/2, maximum(θs), text("θ₀ = $(theta0_deg)°", :center, :top, 10))

    p2 = plot(ts, ls, label="L (length)", xlabel="t", ylabel="L (m)", color=:red, lw=2)
    annotate!(p2, (ts[1] + ts[end])/2, (maximum(ls) + minimum(ls))/2, text("ϕ₀ = $(phi0_deg)°", :center, :center, 10))

    return plot(p1, p2, layout=(2,1))
end
