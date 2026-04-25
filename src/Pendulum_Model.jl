using ModelingToolkit

# ---- Registered symbolic wrappers for MTK ----
# MTK cannot trace through opaque Julia functions (control flow, mod, etc.),
# so we register them. Since we use the explicit λ formula (manual index reduction),
# MTK never needs to symbolically differentiate these — no @register_derivative needed.

@register_symbolic sym_L(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no)
@register_symbolic sym_dL(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no)
@register_symbolic sym_ddL(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no)

sym_L(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no) =
    weightpos(t, l_max, shift, τ, nothing; ϕ₀=ϕ₀, pulse_dist=pulse_dist, g_no=g_no)
sym_dL(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no) =
    weightpos_prim(t, l_max, shift, τ, nothing; ϕ₀=ϕ₀, pulse_dist=pulse_dist, g_no=g_no)
sym_ddL(t, l_max, shift, τ, ϕ₀, pulse_dist, g_no) =
    weightpos_sec(t, l_max, shift, τ, nothing; ϕ₀=ϕ₀, pulse_dist=pulse_dist, g_no=g_no)

"""
    pendulum_xy_model(θ₀, ω₀, tspan, l_max, shift, τ; kwargs...)

Solve the variable-length pendulum in Cartesian (x, y) coordinates using ModelingToolkit.jl.

Uses manually index-reduced equations: the tension λ(t) is computed from the
analytical formula obtained by differentiating the constraint x² + y² = L(t)²
twice (see `TensionIndexReduction.md`). This yields a standard 4-state ODE
(no mass matrix, no DAE), compatible with any explicit solver.

The interface mirrors `solve_pendulum`.
"""
function pendulum_xy_model(θ₀, ω₀, tspan, l_max, shift, τ;
                           ϕ₀=0, pulse_dist=π/2, g_no=2,
                           Q=50, L0=l_max,
                           span=nothing,
                           solver=BS3(), kwargs...)

    @independent_variables t
    @variables x(t) y(t) vx(t) vy(t)
    @variables λ(t) L_val(t) dL_val(t) ddL_val(t)
    @parameters p_l_max p_shift p_τ p_ϕ₀ p_pulse_dist p_g_no p_Q p_L0
    D = Differential(t)

    c_damp = (1 / p_Q) * sqrt(g / p_L0)

    eqs = [
        # Length and its derivatives from registered functions
        L_val ~ sym_L(t, p_l_max, p_shift, p_τ, p_ϕ₀, p_pulse_dist, p_g_no),
        dL_val ~ sym_dL(t, p_l_max, p_shift, p_τ, p_ϕ₀, p_pulse_dist, p_g_no),
        ddL_val ~ sym_ddL(t, p_l_max, p_shift, p_τ, p_ϕ₀, p_pulse_dist, p_g_no),
        # Tension from analytical index reduction (see TensionIndexReduction.md)
        λ ~ (vx^2 + vy^2 - dL_val^2 - L_val*ddL_val - c_damp*L_val*dL_val - g*y) / L_val^2,
        # Equations of motion
        D(x) ~ vx,
        D(y) ~ vy,
        D(vx) ~ -λ * x - c_damp * vx,
        D(vy) ~ -λ * y - c_damp * vy - g,
    ]

    @named sys = ODESystem(eqs, t)
    sys_simp = structural_simplify(sys)

    tspan = Float64.(tspan) # ensure all parameter are Float64

    # Initial conditions: polar → Cartesian
    t0 = tspan[1]
    L_i = weightpos(t0, l_max, shift, τ, span; ϕ₀, pulse_dist, g_no)
    dL_i = weightpos_prim(t0, l_max, shift, τ, span; ϕ₀, pulse_dist, g_no)

    x0  =  L_i * sin(θ₀)
    y0  = -L_i * cos(θ₀)
    vx0 =  dL_i * sin(θ₀) + L_i * ω₀ * cos(θ₀)
    vy0 = -dL_i * cos(θ₀) + L_i * ω₀ * sin(θ₀)

    init = Dict(
        x => x0, 
        y => y0, 
        vx => vx0, 
        vy => vy0,
        p_l_max => Float64(l_max), 
        p_shift => Float64(shift),
        p_τ => Float64(τ), 
        p_ϕ₀ => Float64(ϕ₀),
        p_pulse_dist => Float64(pulse_dist), 
        p_g_no => Float64(g_no),
        p_Q => Float64(Q), 
        p_L0 => Float64(L0),
    )

    prob = ODEProblem(sys_simp, init, tspan)
    return solve(prob, solver; kwargs...)
end

function plot_pendulum_xy_model(; θ₀=-π/4, ω₀=0.0, l_max=3, τ=3.5, shift=0.2,
                                 g_no=1, ϕ₀=0, pulse_dist=π/2, Q=50,
                                 tplot=(0.0, 3*τ), solver=BS3(), kwargs...)
    p_start, p_end = tplot
    tspan = (0.0, Float64(p_end))

    sol = pendulum_xy_model(θ₀, ω₀, tspan, l_max, shift, τ;
                            ϕ₀, pulse_dist, g_no, Q, L0=l_max, solver, kwargs...)

    ts = range(p_start, p_end, length=1000)

    # Reconstruct θ from Cartesian: x = L sinθ, y = -L cosθ  ⟹  θ = atan(x, -y)
    @independent_variables t
    @variables x(t) y(t)
    θs = [atan(sol(t_val, idxs=x), -sol(t_val, idxs=y)) for t_val in ts]
    ls = [weightpos(t_val, l_max, shift, τ; ϕ₀, pulse_dist, g_no) for t_val in ts]

    theta0_deg = round(Int, rad2deg(θ₀))
    phi0_deg   = round(Int, rad2deg(ϕ₀))

    p1 = plot(ts, θs, label="θ (angle)", xlabel="t", ylabel="θ (rad)", lw=2)
    annotate!(p1, (ts[1]+ts[end])/2, maximum(θs),
              text("θ₀ = $(theta0_deg)°", 10))

    p2 = plot(ts, ls, label="L (length)", xlabel="t", ylabel="L (m)", color=:red, lw=2)
    annotate!(p2, (ts[1]+ts[end])/2, (maximum(ls)+minimum(ls))/2,
              text("ϕ₀ = $(phi0_deg)°", 10))

    return plot(p1, p2, layout=(2,1))
end