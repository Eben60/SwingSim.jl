
# ============================================================================
# Triple Pendulum with Angle Control ("Foldable Robo-Gymnast")
#
# Uses Symbolics.jl to derive the Lagrangian equations at module-load time
# and generates optimized numerical functions for M(q,α) and f(q,q̇,α,α̇,α̈).
#
# The ODE is:  M * [q̈₁; q̈₂] = f   (2 free DOFs, α(t) prescribed)
#
# See Triple-Pendulum-Approach.md for the full derivation.
# ============================================================================

# ---- Symbolic derivation (runs once at module load) ----

function _build_gymnast_functions()
    @variables q1_s q2_s dq1_s dq2_s α_s dα_s ddα_s
    @variables l1_s l2_s l3_s m1_s m2_s m3_s g_s b_s

    # Cartesian positions (pivot at origin, y-up, angles from downward vertical)
    x1 = l1_s * sin(q1_s)
    y1 = -l1_s * cos(q1_s)
    x2 = x1 + l2_s * sin(q1_s + q2_s)
    y2 = y1 - l2_s * cos(q1_s + q2_s)
    x3 = x2 + l3_s * sin(q1_s + q2_s + α_s)
    y3 = y2 - l3_s * cos(q1_s + q2_s + α_s)

    xs = [x1, x2, x3]
    ys = [y1, y2, y3]
    ms_sym = [m1_s, m2_s, m3_s]
    all_q = [q1_s, q2_s, α_s]
    dqs = [dq1_s, dq2_s, dα_s]

    # 3×3 generalized inertia matrix: H[a,b] = Σ mᵢ (∂xᵢ/∂qₐ ∂xᵢ/∂qᵦ + ∂yᵢ/∂qₐ ∂yᵢ/∂qᵦ)
    H = [Symbolics.simplify(sum(
            ms_sym[i] * (Symbolics.derivative(xs[i], all_q[a]) * Symbolics.derivative(xs[i], all_q[b]) +
                         Symbolics.derivative(ys[i], all_q[a]) * Symbolics.derivative(ys[i], all_q[b]))
            for i in 1:3))
         for a in 1:3, b in 1:3]

    # Gravity generalized forces: ∂V/∂qₐ = Σ mᵢ g ∂yᵢ/∂qₐ
    Grav = [Symbolics.simplify(sum(ms_sym[i] * g_s * Symbolics.derivative(ys[i], all_q[a]) for i in 1:3))
            for a in 1:3]

    # Coriolis + centrifugal: C_a = Σ_{b,c} Γ_{a,bc} dq_b dq_c
    # where Γ_{a,bc} = ½(∂H_{ab}/∂q_c + ∂H_{ac}/∂q_b - ∂H_{bc}/∂q_a)
    C = Array{Num}(undef, 2)
    for a in 1:2
        val = Num(0)
        for b in 1:3, c in 1:3
            Gamma = (Symbolics.derivative(H[a,b], all_q[c]) +
                     Symbolics.derivative(H[a,c], all_q[b]) -
                     Symbolics.derivative(H[b,c], all_q[a])) / 2
            val += Gamma * dqs[b] * dqs[c]
        end
        C[a] = Symbolics.simplify(val)
    end

    # Mass matrix (2×2, free DOFs only)
    M_sym = H[1:2, 1:2]

    # RHS: M * q̈ = rhs
    #   rhs = -Coriolis - Gravity - H[a,3]*α̈ - damping(q₁ only)
    rhs_sym = [-C[a] - Grav[a] - H[a,3] * ddα_s - (a == 1 ? b_s * dq1_s : Num(0))
               for a in 1:2]

    # Build optimized numerical functions as expressions (safe for precompilation)
    all_vars = [q1_s, q2_s, dq1_s, dq2_s, α_s, dα_s, ddα_s,
                l1_s, l2_s, l3_s, m1_s, m2_s, m3_s, g_s, b_s]

    M_expr = Symbolics.build_function(M_sym, all_vars...)[1]
    rhs_expr = Symbolics.build_function(rhs_sym, all_vars...)[1]

    return M_expr, rhs_expr
end

# Derive symbolic expressions at precompile time
const _gymnast_M_expr, _gymnast_rhs_expr = _build_gymnast_functions()

# Mutable storage for the compiled functions (populated at __init__)
const _gymnast_funcs = Ref{Any}(nothing)

function __init_gymnast__()
    M_f = eval(_gymnast_M_expr)
    rhs_f = eval(_gymnast_rhs_expr)
    _gymnast_funcs[] = (; M_f, rhs_f)
end

# ---- Discontinuity times for tstops ----

"""
    _gymnast_tstops(tspan, shift, τ; ϕ₀=0, pulse_dist=π/2, g_no=2)

Compute all span-boundary times of α(t) within `tspan`.
These are the points where ddα jumps (piecewise-constant second derivative),
so the ODE solver should step exactly to them.
"""
function _gymnast_tstops(tspan, shift, τ; ϕ₀=0, pulse_dist=π/2, g_no=2)
    a = g_no * g
    t_a = sqrt(shift / a)
    t_pulse2 = τ * (pulse_dist / (2π))

    # Discontinuity points within one period [0, τ) in reduced time
    raw = [0.0, t_a, 2t_a, t_pulse2, t_pulse2 + t_a, t_pulse2 + 2t_a]
    boundaries = sort!(unique!(filter(t -> 0 ≤ t < τ, raw)))

    # Map to actual time across all periods covering tspan
    phase_shift = τ * ϕ₀ / (2π)
    t0, tf = tspan
    tstops = Float64[]
    k_min = floor(Int, (t0 + phase_shift) / τ)
    k_max = ceil(Int, (tf + phase_shift) / τ)

    for k in k_min:k_max
        for tb in boundaries
            t_actual = tb - phase_shift + k * τ
            if t0 < t_actual < tf
                push!(tstops, t_actual)
            end
        end
    end

    return sort!(tstops)
end

# ---- ODE function ----

function _gymnast_ode!(du, u, p, t)
    q1, q2, dq1, dq2 = u
    l1, l2, l3, m1, m2, m3, τ, shift, ϕ₀, pulse_dist, g_no, b_coeff = p

    # Prescribed servo angle α(t) and its derivatives
    # Compute span once, pass explicitly to avoid redundant span_no calls
    (; span) = span_no(t, shift, τ; ϕ₀, pulse_dist, g_no)
    α   = weightpos(t, 0.0, shift, τ, span; ϕ₀, pulse_dist, g_no)
    dα  = weightpos_prim(t, 0.0, shift, τ, span; ϕ₀, pulse_dist, g_no)
    ddα = weightpos_sec(t, 0.0, shift, τ, span; ϕ₀, pulse_dist, g_no)

    args = (q1, q2, dq1, dq2, α, dα, ddα, l1, l2, l3, m1, m2, m3, g, b_coeff)

    (; M_f, rhs_f) = _gymnast_funcs[]
    M = M_f(args...)
    rhs = rhs_f(args...)

    # Solve M * [q̈₁; q̈₂] = rhs
    ddq = M \ rhs

    du[1] = dq1      # q̇₁
    du[2] = dq2      # q̇₂
    du[3] = ddq[1]   # q̈₁
    du[4] = ddq[2]   # q̈₂
    return nothing
end

# ---- Solver ----

"""
    solve_pendulum_gymnast(q1₀, q2₀, dq1₀, dq2₀, tspan, l1, l2, l3; kwargs...)

Solve the triple-pendulum (foldable gymnast) model.

# Arguments
- `q1₀, q2₀`: initial angles (rad) for free DOFs
- `dq1₀, dq2₀`: initial angular velocities (rad/s)
- `tspan`: time span tuple
- `l1, l2, l3`: link lengths (m)

# Keyword arguments
- `m1, m2, m3`: point masses (kg), default 1.0
- `τ`: drive period (s), default = natural period of folded configuration
- `shift`: α(t) amplitude, default 0.8π
- `ϕ₀`: phase offset for α(t), default 0
- `pulse_dist`: pulse distribution parameter, default π/2
- `g_no`: acceleration multiplier for α(t) profile
- `Q`: quality factor for damping, default 50
- `solver`: ODE solver, default Tsit5()
"""
function solve_pendulum_gymnast(q1₀, q2₀, dq1₀, dq2₀, tspan, l1, l2, l3;
                                m1=1.0, m2=1.0, m3=1.0,
                                τ=nothing, shift=0.8π,
                                ϕ₀=0, pulse_dist=π/2, g_no=nothing,
                                Q=50.0, solver=Tsit5(), kwargs...)

    # Reference configuration: folded gymnast
    L0 = l1 + l2 - l3
    I0 = m1 * l1^2 + m2 * (l1 + l2)^2 + m3 * L0^2
    ω0 = sqrt(g / L0)

    if isnothing(τ)
        τ = 2π / ω0   # natural period of reference config
    end

    if isnothing(g_no)
        # Minimum acceleration to achieve shift within the period
        a_min = shift * 256 / τ^2
        g_no = round(2 * a_min / g; digits=1)
    end

    # Damping coefficient: b = I₀ ω₀ / Q
    b_coeff = I0 * ω0 / Q

    p = (l1, l2, l3, m1, m2, m3, τ, shift, ϕ₀, pulse_dist, g_no, b_coeff)
    u0 = [q1₀, q2₀, dq1₀, dq2₀]
    tspan = Float64.(tspan)

    # Step exactly to span boundaries where ddα is discontinuous
    stops = _gymnast_tstops(tspan, shift, τ; ϕ₀, pulse_dist, g_no)
    prob = ODEProblem(_gymnast_ode!, u0, tspan, p)
    return solve(prob, solver; tstops=stops, kwargs...)
end

# ---- Plotting ----

"""
    wrap_angle(θ)

Map angle to (-π, π] range.
"""
wrap_angle(θ) = mod(θ + π, 2π) - π

"""
    plot_pendulum_gymnast(; kwargs...)

Test plot for the foldable gymnast model.
Top panel: θ₁ and θ₂ (wrapped to ±π). Bottom panel: θ₃ (servo angle).
"""
function plot_pendulum_gymnast(; q1₀=-π/6, q2₀=0.0, dq1₀=0.0, dq2₀=0.0,
                                 l1=5.0, l2=1.0, l3=0.7,
                                 m1=1.0, m2=1.0, m3=1.0,
                                 Q=50.0, tplot=nothing,
                                 solver=Tsit5(), kwargs...)

    # Compute τ from reference configuration
    L0 = l1 + l2 - l3
    τ = 2π * sqrt(L0 / g)
    shift = 0.8π

    # Minimum acceleration and double it
    a_min = shift * 256 / τ^2
    g_no = round(2 * a_min / g; digits=1)

    if isnothing(tplot)
        tplot = (0.0, 5 * τ)
    end
    p_start, p_end = tplot
    tspan = (0.0, Float64(p_end))

    sol = solve_pendulum_gymnast(q1₀, q2₀, dq1₀, dq2₀, tspan, l1, l2, l3;
                                 m1, m2, m3, τ, shift, Q, g_no, solver, kwargs...)

    ts = range(p_start, p_end, length=1000)

    θ1s = [wrap_angle(sol(t)[1]) for t in ts]
    θ2s = [wrap_angle(sol(t)[2]) for t in ts]

    # α(t) from weightpos, which is the servo angle θ₃
    αs = [weightpos(t, 0.0, shift, τ; g_no) for t in ts]

    q1_deg = round(Int, rad2deg(q1₀))

    p1 = plot(ts, θ1s, label="θ₁", xlabel="t (s)", ylabel="angle (rad)", lw=2)
    plot!(p1, ts, θ2s, label="θ₂", lw=2)
    annotate!(p1, (ts[1] + ts[end]) / 2, maximum(θ1s),
              text("q₁₀ = $(q1_deg)°", 10))

    p2 = plot(ts, αs, label="θ₃", xlabel="t (s)", ylabel="angle (rad)",
              color=:purple, lw=2)

    return plot(p1, p2, layout=(2, 1))
end
