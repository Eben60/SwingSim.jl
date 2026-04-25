
function _weightpos(t, l_max, shift, τ, span; pulse_dist = π/2, g_no)
    a = g_no * g
    t_a = sqrt(shift / a)
    t_pulse2 = τ * (pulse_dist / (2π))

    dl = if span == 1
        0.5 * a * t^2
    elseif span == 2
        0.5 * shift
    elseif span == 3
        dt = t - t_a
        shift/2 + a * dt * (t_a - dt/2)
    elseif span == 4
        shift
    elseif span == 5
        dt = t - t_pulse2
        shift - 0.5 * a * dt^2
    elseif span == 6
        0.5 * shift
    elseif span == 7
        dt = t - (t_pulse2 + t_a)
        shift/2 - a * dt * (t_a - dt/2)
    else
        0.0
    end
    return l_max - dl
end

τ_so(l) = 2π * sqrt(l / g)

# function weightpos(t, l_max, shift, τ; ϕ₀=0, pulse_dist = π/2, g_no = 2)
#     (; span) = span_no(t, shift, τ; ϕ₀, pulse_dist, g_no)
#     return weightpos(t, l_max, shift, τ, span; ϕ₀, pulse_dist, g_no)
# end

function weightpos(t, l_max, shift, τ, span=nothing; ϕ₀=0, pulse_dist = π/2, g_no = 2)
    t_red = t_reduced(t, τ, ϕ₀)
    if isnothing(span) 
        (; span) = span_no(t, shift, τ; ϕ₀, pulse_dist, g_no)
    end
    return _weightpos(t_red, l_max, shift, τ, span; pulse_dist, g_no)
end

function _weightpos_prim(t, l_max, shift, τ, span; pulse_dist = π/2, g_no)
    a = g_no * g
    t_a = sqrt(shift / a)
    t_pulse2 = τ * (pulse_dist / (2π))

    ddl = if span == 1
        a * t
    elseif span == 2
        0.0
    elseif span == 3
        dt = t - t_a
        a * (t_a - dt)
    elseif span == 4
        0.0
    elseif span == 5
        dt = t - t_pulse2
        -a * dt
    elseif span == 6
        0.0
    elseif span == 7
        dt = t - (t_pulse2 + t_a)
        -a * (t_a - dt)
    else
        0.0
    end
    return -ddl
end

function weightpos_prim(t, l_max, shift, τ, span=nothing; ϕ₀=0, pulse_dist = π/2, g_no = 2)
    t_red = t_reduced(t, τ, ϕ₀)
    if isnothing(span) 
        (; span) = span_no(t, shift, τ; ϕ₀, pulse_dist, g_no)
    end
    return _weightpos_prim(t_red, l_max, shift, τ, span; pulse_dist, g_no)
end


t_reduced(t, τ, ϕ₀=0) = mod(t + τ*ϕ₀/(2*π), τ)

span_no(t, shift, τ; ϕ₀=0, pulse_dist = π/2, g_no = 2) =
    _span_no(t_reduced(t, τ, ϕ₀) , shift, τ, 10*eps(float(t)) ;  pulse_dist, g_no )


function _span_no(t, shift, τ, tol;  pulse_dist, g_no )
    a = g_no * g
    # Maximum achievable shift within τ/8 total time (t_a = τ/16) is a * (τ/16)^2
    shift > a * (τ/16)^2 && error("Shift $shift cannot be achieved within time at acceleration $a")

    t_a = sqrt(shift / a)
    t_pulse2 = τ * (pulse_dist / (2π))

    t1 = t_a
    t2 = t1
    t3 = t2 + t_a
    t4 = t_pulse2
    t5 = t4 + t_a
    t6 = t5
    t7 = t6 + t_a

    if t < t1 - tol
        return (;span=1, next_starting=t1)
    elseif t < t2 - tol
        return (;span=2, next_starting=t2)
    elseif t < t3 - tol
        return (;span=3, next_starting=t3)
    elseif t < t4 - tol
        return (;span=4, next_starting=t4)
    elseif t < t5 - tol
        return (;span=5, next_starting=t5)
    elseif t < t6 - tol
        return (;span=6, next_starting=t6)
    elseif t < t7 - tol
        return (;span=7, next_starting=t7)
    else
        return (;span=8, next_starting=τ)
    end
end


function weightspos_testplot(; l_max=3, τ=3.5, shift=0.2, g_no=1, ϕ₀=0, pulse_dist=π/2)
    ts = range(0, 2τ, length=1000)
    ys = Float64[weightpos(t, l_max, shift, τ; ϕ₀, pulse_dist, g_no) for t in ts]
    return plot(ts, ys)
end

# ---- Second derivative of length (L̈) ----

function _weightpos_sec(t, l_max, shift, τ, span; pulse_dist = π/2, g_no)
    a = g_no * g

    ddl = if span == 1
        a
    elseif span == 2
        0.0
    elseif span == 3
        -a
    elseif span == 4
        0.0
    elseif span == 5
        -a
    elseif span == 6
        0.0
    elseif span == 7
        a
    else
        0.0
    end
    return -ddl
end

function weightpos_sec(t, l_max, shift, τ, span=nothing; ϕ₀=0, pulse_dist = π/2, g_no = 2)
    t_red = t_reduced(t, τ, ϕ₀)
    if isnothing(span)
        (; span) = span_no(t, shift, τ; ϕ₀, pulse_dist, g_no)
    end
    return _weightpos_sec(t_red, l_max, shift, τ, span; pulse_dist, g_no)
end