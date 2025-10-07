##############################
# Imports
##############################
using HomotopyContinuation
using Random
using LinearAlgebra
using ImplicitPlots, Plots  # Plots backend for ImplicitPlots
using IntervalSets             # a..b interval syntax for implicit_plot

# -------------------------------------------------------------------
# Global symbols needed by J1.jl and J2.jl (must be visible to `Main`)
# -------------------------------------------------------------------
using HomotopyContinuation

@var A11 A12 A21 A22 B111 B112 B122 B211 B212 B222 r1 r2

# Load the boundary polynomials once (now they can “see” the symbols)
const J1_POLY = include(joinpath(@__DIR__, "J1.jl"))  # zero/complex boundary
const J2_POLY = include(joinpath(@__DIR__, "J2.jl"))  # bifurcation boundary


##############################
# Existing content (kept)
##############################

############SAMPLE PARAMETERS#############

# Constraint matrix builder
"""
    build_C(ℓ, m, n)

Construct the 2n×D matrix C enforcing:
 1) For each i=1:n, sum over A[i, j] plus r[i] == 0
 2) For each i=1:n, sum over B[:,:,i] plus r[i] == 0
Where D = n + m*n + ℓ*m*n.
"""
function build_C(ℓ::Int, m::Int, n::Int)
    D = n + m*n + ℓ*m*n
    C = zeros(Float64, 2n, D)

    # 1) sum_j A[i,j] + r[i] == 0  (row sums of A)
    for i in 1:n
        C[i, i] = 1.0
        for j in 1:n
            col = n + (j-1)*m + i
            C[i, col] = 1.0
        end
    end

    # 2) sum_{j,k} B[j,k,i] + r[i] == 0
    offsetB = n + m*n
    for i in 1:n
        row = n + i
        C[row, i] = 1.0
        for k in 1:m, j in 1:ℓ
            lin = j + (k-1)*ℓ + (i-1)*ℓ*m
            C[row, offsetB + lin] = 1.0
        end
    end

    return C
end

# Modified hit-and-run sampler using a truncated Normal(0,1) step
"""
    hit_and_run_sample(C, L, U, x0, ℓ, m, n)

Perform one hit-and-run move sampling λ from a truncated Normal(0,1) on the feasible chord.
Returns (r, A, B)
"""
function hit_and_run_sample(
    C::AbstractMatrix{<:Float64},
    L::AbstractVector{<:Float64},
    U::AbstractVector{<:Float64},
    x0::AbstractVector{<:Float64},
    ℓ::Int, m::Int, n::Int;
    max_tries::Int = 100
)
    # Precompute nullspace basis
    F     = svd(C; full=true)
    tol   = maximum(size(C)) * eps(eltype(C)) * maximum(F.S)
    rankC = count(>(tol), F.S)
    Q     = F.V[:, rankC+1:end]
    D     = length(x0)

    for attempt in 1:max_tries
        # random direction in ker(C)
        u = Q * randn(size(Q,2))
        nu = norm(u)
        if nu < 1e-12
            continue
        end
        u ./= nu

        # find feasible interval
        λ_min, λ_max = -Inf, Inf
        for j in 1:D
            if abs(u[j]) > 1e-12
                a = (L[j] - x0[j]) / u[j]
                b = (U[j] - x0[j]) / u[j]
                lo, hi = min(a,b), max(a,b)
                λ_min = max(λ_min, lo)
                λ_max = min(λ_max, hi)
            end
        end
        if λ_max < λ_min || !isfinite(λ_min) || !isfinite(λ_max)
            continue
        end

        # sample λ from truncated Normal(0,1)
        λ = rand() * (λ_max - λ_min) + λ_min
        x = x0 .+ λ .* u

        # split back into (r, A, B)
        r =     x[1                 : n]
        A = reshape(x[n+1             : n + m*n], m, n)
        B = reshape(x[n + m*n + 1     : end],       ℓ, m, n)

        return r, A, B
    end
    error("Failed to find a valid hit-and-run step after $max_tries tries.")
end

"""
    generate_parameter_set(ℓ::Int, m::Int, n::Int;
                           lower::Float64 = -10.0,
                           upper::Float64 = 10.0,
                           seed::Union{Nothing,Int}=nothing)

Run one hit-and-run sample and return the tuple `(r, A, B)`.
No files are written.
"""
function generate_parameter_set(ℓ::Int, m::Int, n::Int;
                                lower::Float64 = -10.0,
                                upper::Float64 = 10.0,
                                seed::Union{Nothing,Int}=nothing)
    if seed !== nothing
        Random.seed!(seed)
    end

    D = n + m*n + ℓ*m*n
    C = build_C(ℓ, m, n)
    L = fill(lower, D)
    U = fill(upper, D)
    x0 = vcat(zeros(n), zeros(m*n), zeros(ℓ*m*n))

    r, A, B = hit_and_run_sample(C, L, U, x0, ℓ, m, n)
    return r, A, B
end


###########FORM SYSTEM#############
"""
    build_num_eqs(r, A, B, x) -> Vector{Expression}

Given sampled parameters:
  • r :: AbstractVector{<:Real}         (length n)
  • A :: AbstractMatrix{<:Real}         (n × n)
  • B :: AbstractArray{<:Real,3}        (n × n × n)

and variables:
  • x :: Vector{Variable}               (length n)

return the vector of numerical polynomials `num_eqs` with:
    f_i(x) = r[i] + ∑_j A[i,j] * x[j] + ∑_{j,k} B[i,j,k] * x[j] * x[k]
"""
function build_num_eqs(r::AbstractVector, A::AbstractMatrix, B::AbstractArray{<:Real,3}, x::Vector{Variable})
    n = length(x)
    @assert length(r) == n "length(r) must equal length(x)"
    @assert size(A) == (n, n) "A must be n×n"
    @assert size(B) == (n, n, n) "B must be n×n×n (i,j,k)"

    eqs = Vector{Expression}(undef, n)
    for i in 1:n
        lin  = sum(A[i, j] * x[j] for j in 1:n)
        quad = sum(B[j, k, i] * x[j] * x[k] for j in 1:n, k in 1:n)
        eqs[i] = r[i] + lin + quad
    end
    return eqs
end

"""
    get_parametrized_system(num_eqs::Vector{Expression},
                            x::Vector{Variable},
                            Δr::Vector{Variable},
                            α::Variable)

Given numerical polynomials `num_eqs` in variables `x`, return a parameterized `System` where:
  • constant term r0ᵢ → r0ᵢ + Δrᵢ
  • linear terms (order 1) → (1 - α) * coeff
  • quadratic terms (order 2) → α * coeff
Assumes degree ≤ 2. Parameters are `vcat(Δr, α)`.
"""
function get_parametrized_system(num_eqs::Vector{Expression},
                                 x::Vector{Variable},
                                 Δr::Vector{Variable},
                                 α::Variable)
    @assert length(num_eqs) == length(Δr) "Need one Δr per equation/species."

    out = Vector{Expression}(undef, length(num_eqs))

    for i in eachindex(num_eqs)
        f = num_eqs[i]
        exps, coeffs = exponents_coefficients(f, x)
        orders = vec(sum(exps; dims = 1))  # values in {0,1,2}

        terms = Vector{Expression}(undef, length(coeffs))
        @inbounds for k in eachindex(coeffs)
            # build the monomial from its exponent vector (matches coeffs[k])
            mon = prod(x[j]^exps[j,k] for j in 1:length(x))

            # reweight coefficient by order
            c = coeffs[k]
            if orders[k] == 0
                c = c + Δr[i]         # r0ᵢ → r0ᵢ + Δrᵢ
            elseif orders[k] == 1
                c = (1 - α) * c       # linear
            else # orders[k] == 2
                c = α * c             # quadratic
            end

            terms[k] = c * mon
        end

        out[i] = sum(terms)
    end

    return System(out; variables = x, parameters = vcat(Δr, α))
end

#####################BUILD NEGATIVE AND COMPLEX BOUNDARY EQUATIONS######################

function build_boundary_system(r0::AbstractVector, A::AbstractMatrix, B::AbstractArray{<:Real,3};
                               dir::AbstractString = @__DIR__)
    @assert length(r0) == 2
    @assert size(A) == (2,2)
    @assert size(B) == (2,2,2)

    @var a Δr1 Δr2

    # Start from globally-loaded polynomials
    g1 = J1_POLY
    g2 = J2_POLY

    # Symbols & numeric values
    Asyms = (A11, A12, A21, A22)
    Bsyms = (B111, B112, B122, B211, B212, B222)
    Avals = (A[1,1], A[1,2], A[2,1], A[2,2])
    Bvals = (B[1,1,1], B[1,1,2], B[1,2,2], B[2,1,1], B[2,1,2], B[2,2,2])

    # Scale A,B: A ↦ (1-α)A, B ↦ αB  (varargs: splat pairs)
    g1 = subs(g1, (Asyms .=> ((1 - a) .* Asyms))..., (Bsyms .=> (a .* Bsyms))...)
    g2 = subs(g2, (Asyms .=> ((1 - a) .* Asyms))..., (Bsyms .=> (a .* Bsyms))...)

    # r = r0 + Δr  (use varargs pairs, not tuple→tuple)
    g1 = subs(g1, r1 => r0[1] + Δr1, r2 => r0[2] + Δr2)
    g2 = subs(g2, r1 => r0[1] + Δr1, r2 => r0[2] + Δr2)

    # Now plug numeric A,B (again, varargs of pairs)
    g1Δ = subs(g1, (Asyms .=> Avals)..., (Bsyms .=> Bvals)...)
    g2Δ = subs(g2, (Asyms .=> Avals)..., (Bsyms .=> Bvals)...)

    fΔ = System([g1Δ; g2Δ]; variables=[Δr1, Δr2], parameters=[a])
    return fΔ, g1Δ, g2Δ
end


##############################
# NEW: equilibria & plotting
##############################

"""
    compute_equilibria(r0, A, B; α, real_only=true, tol=1e-8)

Solve f(x) = 0 for Δr = 0 and a given α.
Returns (sols, result), where `sols` are real solutions (if `real_only=true`).
"""
function compute_equilibria(r0::AbstractVector, A::AbstractMatrix, B::AbstractArray{<:Real,3};
                            α::Real, real_only::Bool=true, tol::Real=1e-8)
    n = length(r0)
    @var x[1:n] Δr[1:n] a

    # Build base polynomials and parameterize with Δr, α
    num_eqs = build_num_eqs(r0, A, B, x)
    S = get_parametrized_system(num_eqs, x, Δr, a)

    # Evaluate Δr = 0 and a = α (do substitutions sequentially)
    eqs = [subs(eq, Δr => zeros(n)) for eq in S.expressions]
    eqs = [subs(eq, a => α) for eq in eqs]

    Sα = System(eqs; variables = x)
    R = solve(Sα)
    sols = solutions(R)

    if real_only
        sols = [real.(z) for z in sols if maximum(abs.(imag.(z))) ≤ tol]
    end
    return sols, R
end


function plot_boundaries_implicit(r0, A, B; α::Real=0.3,
                                  xlims::Tuple{<:Real,<:Real}=(-2,2),
                                  ylims::Tuple{<:Real,<:Real}=(-2,2))
    fΔ, g1Δ, g2Δ = build_boundary_system(r0, A, B)
    a_sym = parameters(fΔ)[1]
    g1α = subs(g1Δ, a_sym => α)
    g2α = subs(g2Δ, a_sym => α)

    p = implicit_plot(g1α; xlims=xlims, ylims=ylims, linewidth=2, label="zero (J1)", linecolor = :indianred,
                      xlabel="Δr₁", ylabel="Δr₂", title="α=$(round(α,digits=3))")
    implicit_plot!(p, g2α; xlims=xlims, ylims=ylims, linewidth=2, label="bif (J2)", linecolor = :steelblue)
    savefig(p, "boundary_plot.pdf")
end

using ImplicitPlots, Plots

"""
    save_boundaries_gif_implicit(r0, A, B;
                                 alphas = 0.1:0.1:0.9,
                                 xlims = (-2.0, 2.0),
                                 ylims = (-2.0, 2.0),
                                 outfile = "boundary_anim.gif",
                                 fps = 2)

Animate the J1 (zero) and J2 (bifurcation) boundaries in (Δr₁, Δr₂) as α varies.
Writes a GIF to `outfile`. Returns nothing.
"""
function save_boundaries_gif_implicit(r0, A, B;
                                      alphas = 0.1:0.1:0.9,
                                      xlims = (-2.0, 2.0),
                                      ylims = (-2.0, 2.0),
                                      outfile = "boundary_anim.gif",
                                      fps::Integer = 2)
    # Build boundary expressions
    fΔ, g1Δ, g2Δ = build_boundary_system(r0, A, B)
    a_sym = parameters(fΔ)[1]

    # Optional: plot defaults
    default(size=(700, 540), legend=:topright, framestyle=:box)

    anim = @animate for ai in alphas
        println("Rendering α = $ai ...")
        g1a = subs(g1Δ, a_sym => ai)
        g2a = subs(g2Δ, a_sym => ai)

        implicit_plot(g1a;
            xlims=xlims, ylims=ylims,
            linewidth=2, linecolor=:indianred, label="zero boundary",
            xlabel="Δr₁", ylabel="Δr₂",
            title="α = $(round(ai, digits=3))"
        )
        implicit_plot!(g2a;
            xlims=xlims, ylims=ylims,
            linewidth=2, linecolor=:steelblue, label="bifurcation"
        )
    end

    gif(anim, outfile; fps=fps)
    return nothing
end

"""
given start and end parameters, and a value of t, return the corresponding convex sum parameters
"""
function get_parameters_at_t(t::Float64, initial_parameters::Vector{Float64}, target_parameters::Vector{Float64})
    t * initial_parameters + (1 - t) * target_parameters
end


function decide_boundary_type(
    syst::System,
    tracker::Tracker,
    result::TrackerResult,
    tol::Float64,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64}
)
    return_code = result.return_code
    x = result.solution

    if return_code == :success
        # Tracking was successful and solution is positive
        return :boundary

    elseif return_code == :tracking
        # Tracking was interrupted manually by trackpositive!
        if any(abs.(real.(x)) .< tol)
            #println("Solution at negative boundary: ", x)
            return :negative
        else
            #println("Solution became negative but above tolerance: ", x)
            return :nonconverged
        end

    elseif return_code == :terminated_step_size_too_small || return_code == :terminated_max_steps
        # Tracking didn't succeed nor found negative solutions only option left is complex bifurcation.
           return :complex
    else 
        # I don't know what this could be: print
        #println("Unkown boundary type: ", return_code)
        return :unknown
    end
end

# --- Tracking with positivity constraint ---

"""
    trackpositive!(tracker::Tracker, x::AbstractVector, [t1=1.0, t0=0.0]; kwargs...) -> (t_before, x_before)

Track a solution path, stop when any component becomes negative.
"""
function trackpositive!(
    tracker::Tracker,
    x::AbstractVector,
    t1 = 1.0,
    t0 = 0.0;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    extended_precision::Bool = false,
    τ::Float64 = Inf,
    keep_steps::Bool = false,
    max_initial_step_size::Float64 = Inf,
    debug::Bool = false,
)
    init!(
        tracker,
        x,
        t1,
        t0;
        ω = ω,
        μ = μ,
        extended_precision = extended_precision,
        τ = τ,
        keep_steps = keep_steps,
        max_initial_step_size = max_initial_step_size,
    )
    #record initial t and equilibrium
    tbefore = t1
    xbefore = x
    while is_tracking(tracker.state.code)
        #record quantities before and after step
        tbefore = tracker.state.t
        xbefore = copy(tracker.state.x)
        step!(tracker, debug)
        #after a tracking step, check if any of the components are negative.
        if any(real(tracker.state.x) .< 0)
            #return the component before negative solution was found
            return real(tbefore), real(xbefore)
        end
    end
    return real(tbefore), real(xbefore)
end

# --- Higher level tracking monitor ---

"""
    track_and_monitor!(syst, initialsol, initial_parameters, target_parameters, tol) -> (t_before, x_before, t_after, x_after, flag)

Track once from initial to target parameters, detect feasibility loss.
Returns last positive solution, solution after tracking, and boundary flag.
"""
function track_and_monitor!(
    syst::System,
    initialsol::AbstractVector,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64},
    tol::Float64,
    max_step_ratio::Float64 = 1.
)   
    par_dist = norm(initial_parameters - target_parameters)
    #the maximum step size should be 10 times smaller than the parameter euclidean distance
    max_step_size = par_dist * max_step_ratio 

    tracker = Tracker(CoefficientHomotopy(syst; start_coefficients = initial_parameters,
                                           target_coefficients = target_parameters);
                                           options = TrackerOptions(max_steps = Int(1e6))) #increased number of steps
    t_before, x_before = trackpositive!(tracker, initialsol, 1.0, 0.0, max_initial_step_size = max_step_size) #adjusted max_initial_step_size
    result = TrackerResult(tracker.homotopy, tracker.state)
    t_after = real(result.t)
    x_after = result.solution
    #store some tracking output
    flag = decide_boundary_type(syst, tracker, result, tol, initial_parameters, target_parameters)
    if flag == :complex
        # println("AT COMPLEX BOUNDARY")
        # println("t_before: ", t_after)
    end
    return t_before, x_before, t_after, x_after, flag
end


# --- Output formatter ---
"""
    process_output_boundary(pars_crit, xstar_crit, flag)

Formats output after finding a critical point.

Returns a NamedTuple with:
- flag: Symbol (:success, :negative, :complex, or :nonconverged)
- pars_crit: Vector of critical parameters
- xstar_crit: Vector of equilibrium at boundary
"""
function process_output_boundary(
    pars_crit::Vector{Float64},
    xstar_crit::Vector{ComplexF64},
    flag::Symbol
)
    return (
        pars_crit = pars_crit,
        xstar_crit = xstar_crit,
        flag = flag
    )
end

# --- Recursive critical parameter finder ---

function findparscrit(
    syst::System,
    initialsol::AbstractVector,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64},
    tol::Float64 = 1e-9,
    rec_level::Int = 1
)
    #println("Recursion level: ", rec_level)
    #println("Initial parameters: ", initial_parameters)
    #println("Target parameters: ", target_parameters)
    t_before, x_before, t_after, x_after, flag = track_and_monitor!(
        syst, initialsol, initial_parameters, target_parameters, tol)
    #get initial and end parameters after tracking
    initpars = get_parameters_at_t(t_before, initial_parameters, target_parameters)
    endpars = get_parameters_at_t(t_after, initial_parameters, target_parameters)
    #println("Flag: ", flag)
    if flag == :negative || flag == :complex || flag == :boundary
        # No boundary was crossed, max perturbation reached
        output = process_output_boundary(endpars, x_after, flag)
        #println("")
        #println("Finished at: ", flag, " equilibrium: ", output.xstar_crit)
        return output
    elseif rec_level > 5
        #maximum recursion level reached, return result
        output = process_output_boundary(target_parameters, x_after, flag)
        return output
    else
        return findparscrit(syst, x_before, initpars, endpars, tol, rec_level + 1)
    end
end

using LinearAlgebra: norm
using Random
using ImplicitPlots, Plots  # only needed if you also plot later

"""
    grid_perts2D_perimeter(; halfwidth=2.0, n_per_side=101)

Return boundary points of a square centered at 0 with side length 2*halfwidth.
Points are ordered clockwise starting at (-halfwidth, +halfwidth).
Corners are included once, so the length is 4*n_per_side - 4.

Example (size 4 square): `grid_perts2D_perimeter(halfwidth=2.0, n_per_side=101)`
"""
function grid_perts2D_perimeter(; halfwidth::Real=2.0, n_per_side::Int=101)
    @assert n_per_side ≥ 2 "n_per_side must be at least 2"
    W = float(halfwidth)

    xs = collect(range(-W,  W;  length=n_per_side))
    ys = collect(range( W, -W;  length=n_per_side))  # top→bottom for clockwise ordering

    perim = NTuple{2,Float64}[]

    # Top edge:    (-W, W) → ( W, W)      (include both corners here)
    append!(perim, ((x,  W) for x in xs))

    # Right edge:  ( W, W) → ( W,-W)      (skip top-right corner to avoid dup)
    append!(perim, (( W, y) for y in ys[2:end]))

    # Bottom edge: ( W,-W) → (-W,-W)      (skip bottom-right corner)
    append!(perim, ((x, -W) for x in Iterators.reverse(xs[1:end-1])))

    # Left edge:   (-W,-W) → (-W, W)      (skip both corners)
    append!(perim, ((-W, y) for y in ys[2:end-1]))

    return perim
end


# ---------- build param system f(x; Δr, a) once ----------
"Build the GLV–HOI system f(x; Δr, α) from (r0, A, B). Assumes n=2."
function build_param_system(r0::AbstractVector, A::AbstractMatrix, B::AbstractArray{<:Real,3})
    n = length(r0)
    @assert n == 2 "This helper is 2D-only."
    @var x[1:n] Δr[1:n] a
    num_eqs = build_num_eqs(r0, A, B, x)
    S = get_parametrized_system(num_eqs, x, Δr, a)   # params: (Δr1, Δr2, a)
    return S, x, Δr, a
end

using ImplicitPlots, Plots
using LinearAlgebra: norm

using CairoMakie

"""
    implicit_contours_makie(g1α, g2α, Δr_vars;
                            xlims=(-2,2), ylims=(-2,2), N=401,
                            colors=(:indianred, :steelblue),
                            outfile=nothing)

Makie equivalent of ImplicitPlots:
evaluate g1α, g2α on a grid in (Δr₁, Δr₂) and draw their 0-contours.

- g1α, g2α :: Expression           # after substituting α => value
- Δr_vars  :: Vector{Variable}     # [Δr1, Δr2] in that order
- outfile  :: String or `nothing`  # if set, figure is saved
"""
function implicit_contours_makie(g1α, g2α, Δr_vars;
                                 xlims=(-2,2), ylims=(-2,2), N=401,
                                 colors=(:indianred, :steelblue),
                                 outfile=nothing)

    xs = range(Float64(xlims[1]), Float64(xlims[2]); length=N)
    ys = range(Float64(ylims[1]), Float64(ylims[2]); length=N)

    # evaluate on a grid
    Z1 = [evaluate(g1α, Δr_vars[1] => x, Δr_vars[2] => y) for x in xs, y in ys]
    Z2 = [evaluate(g2α, Δr_vars[1] => x, Δr_vars[2] => y) for x in xs, y in ys]

    fig = CairoMakie.Figure(resolution=(860, 620))
    ax  = CairoMakie.Axis(fig[1,1], xlabel="Δr₁", ylabel="Δr₂")

    # draw 0-contours (transpose grid to match Makie’s orientation)
    c1 = CairoMakie.contour!(ax, xs, ys, Z1; levels=[0], color=colors[1], linewidth=2)
    c2 = CairoMakie.contour!(ax, xs, ys, Z2; levels=[0], color=colors[2], linewidth=2)
    axislegend(ax, [c1, c2], ["zero boundary", "bifurcation"])

    if outfile !== nothing
        CairoMakie.save(outfile, fig)
    end
    return nothing
end

using HomotopyContinuation
using ImplicitPlots, Plots

using HomotopyContinuation
using ImplicitPlots, Plots

"""
    plot_custom_boundaries_implicit(r0; α=0.3, xlims=(-2,2), ylims=(-2,2),
                                    outfile="custom_boundaries.pdf")

Plot the zero-level sets of the two *explicit* boundary polynomials g1 and g2
(in variables r1, r2, a), mapped to the (Δr₁, Δr₂) plane via r = r0 + Δr and a = α.
- g1 in indianred, g2 in steelblue.
"""
function plot_custom_boundaries_implicit(r0::AbstractVector;
        α::Real=0.3, xlims::Tuple{<:Real,<:Real}=(-2,2),
        ylims::Tuple{<:Real,<:Real}=(-2,2),
        outfile::AbstractString="custom_boundaries.pdf")

    @assert length(r0) == 2 "r0 must be length 2."

    @var r1 r2 a Δr1 Δr2

    # ---- g1 (your previous explicit polynomial) ----
    g1 = 198450*r1^3*a - 1045800*r1^2*r2*a + 1088000*r1*r2^2*a - 320000*r2^3*a -
         7945*r1^2*a^2 + 82628*r1*r2*a^2 - 203392*r2^2*a^2 +
         15890*r1^2*a - 165256*r1*r2*a + 406784*r2^2*a -
         7945*r1^2 + 82628*r1*r2 - 203392*r2^2

    # ---- g2 (your new explicit polynomial) ----
    g2 = 47157273856000000*r1^4*a^4 - 867203893248000000*r1^3*r2*a^4 +
         4497653826048000000*r1^2*r2^2*a^4 - 4696415889408000000*r1*r2^3*a^4 +
         1383054679296000000*r2^4*a^4 + 305662737826880000*r1^3*a^5 -
         926151988939680000*r1^2*r2*a^5 - 2383509908520000000*r1*r2^2*a^5 +
         2300931414228960000*r2^3*a^5 - 356324408034277425*r1^2*a^6 +
         1510874376967760400*r1*r2*a^6 - 339832740214224000*r2^2*a^6 +
         79096139359061580*r1*a^7 - 304298354628175824*r2*a^7 +
         297084481574832*a^8 - 611325475653760000*r1^3*a^4 +
         1852303977879360000*r1^2*r2*a^4 + 4767019817040000000*r1*r2^2*a^4 -
         4601862828457920000*r2^3*a^4 + 1425297632137109700*r1^2*a^5 -
         6043497507871041600*r1*r2*a^5 + 1359330960856896000*r2^2*a^5 -
         474576836154369480*r1*a^6 + 1825790127769054944*r2*a^6 -
         2376675852598656*a^7 + 305662737826880000*r1^3*a^3 -
         926151988939680000*r1^2*r2*a^3 - 2383509908520000000*r1*r2^2*a^3 +
         2300931414228960000*r2^3*a^3 - 2137946448205664550*r1^2*a^4 +
         9065246261806562400*r1*r2*a^4 - 2038996441285344000*r2^2*a^4 +
         1186442090385923700*r1*a^5 - 4564475319422637360*r2*a^5 +
         8318365484095296*a^6 + 1425297632137109700*r1^2*a^3 -
         6043497507871041600*r1*r2*a^3 + 1359330960856896000*r2^2*a^3 -
         1581922787181231600*r1*a^4 + 6085967092563516480*r2*a^4 -
         16636730968190592*a^5 - 356324408034277425*r1^2*a^2 +
         1510874376967760400*r1*r2*a^2 - 339832740214224000*r2^2*a^2 +
         1186442090385923700*r1*a^3 - 4564475319422637360*r2*a^3 +
         20795913710238240*a^4 - 474576836154369480*r1*a^2 +
         1825790127769054944*r2*a^2 - 16636730968190592*a^3 +
         79096139359061580*r1*a - 304298354628175824*r2*a +
         8318365484095296*a^2 - 2376675852598656*a + 297084481574832

    # Map to Δr-coordinates: r = r0 + Δr, then fix α
    g1Δα = subs(subs(g1, r1 => r0[1] + Δr1, r2 => r0[2] + Δr2), a => α)
    g2Δα = subs(subs(g2, r1 => r0[1] + Δr1, r2 => r0[2] + Δr2), a => α)

    # Plot both implicit curves with Plots/ImplicitPlots
    default(size=(860, 620), legend=:topright, framestyle=:box)
    p = implicit_plot(g1Δα; xlims=xlims, ylims=ylims, linewidth=2,
                      linecolor=:indianred, label="g1 (zero boundary)",
                      xlabel="Δr₁", ylabel="Δr₂",
                      title="α = $(round(α, digits=3))")
    implicit_plot!(p, g2Δα; xlims=xlims, ylims=ylims, linewidth=2,
                   linecolor=:steelblue, label="g2 (blue)")

    mkpath(dirname(outfile))
    Plots.savefig(p, outfile)
    return nothing
end



"""
    save_boundary_overlay_perimeter(r0, A, B;
        α=0.3, halfwidth=2.0, n_per_side=101,
        xstar_target=(1.0,1.0),
        outfile="boundary_perimeter_points.pdf")

Plots the implicit boundaries (J1/J2) in (Δr₁, Δr₂) for fixed α, then overlays
numeric boundary points computed along rays from (0,0) to each point on the
perimeter of a square of side 2*halfwidth centered at 0 (regular sampling with
n_per_side points per edge). Points are colored by decide_boundary_type-style
flags (:boundary, :negative, :nonconverged, :complex, :unknown).

Nothing is saved if the implicit plot fails.
"""
function save_boundary_overlay_perimeter(r0::AbstractVector, A::AbstractMatrix, B::AbstractArray{<:Real,3};
                                         α::Real=0.3,
                                         halfwidth::Real=2.0,
                                         n_per_side::Int=101,
                                         xstar_target::Tuple{<:Real,<:Real}=(1.0,1.0),
                                         outfile::AbstractString="boundary_perimeter_points.pdf")

    @assert length(r0) == 2 "This function assumes a 2D system."

    # --- 1) Implicit boundaries from J1/J2 ---
    fΔ, g1Δ, g2Δ = build_boundary_system(r0, A, B)   # variables: [Δr1,Δr2], parameter: [a]
    Δr_vars = variables(fΔ)
    a_sym = parameters(fΔ)[1]
    g1α = subs(g1Δ, a_sym => α)
    g2α = subs(g2Δ, a_sym => α)

    x_lims = (-float(halfwidth), float(halfwidth))
    y_lims = x_lims

    implicit_contours_makie(g1α, g2α, Δr_vars; xlims=(-5,5), ylims=(-5,5), N=401, outfile="boundary_plot.pdf")

    # Build system and get real equilibria at Δr=0, α
    S, x, Δr_syms, a_sym = build_param_system(r0, A, B)
    sols, _ = compute_equilibria(r0, A, B; α=α, real_only=true)
    isempty(sols) && error("No real equilibrium found at Δr=0, α=$α.")
    # Keep only positive equilibria (componentwise > pos_tol)
    pos_sols = [s for s in sols if all(s .> 0.0)]

    # Pick x⋆ as the positive equilibrium closest to xstar_target
    target = collect(xstar_target)
    best_idx = argmin(norm(s .- target) for s in pos_sols)
    xstar = pos_sols[best_idx]

    # ... (rest of your function unchanged)
    p0 = [0.0, 0.0, α]
    perts = grid_perts2D_perimeter(; halfwidth=(x_lims[2]-x_lims[1])/2, n_per_side=n_per_side)  # or use your chosen n_per_side

    points = Float64[]
    flags  = Symbol[]
    sizehint!(points, 2*length(perts))
    sizehint!(flags, length(perts))

    for (δ1, δ2) in perts
        p1 = [δ1, δ2, α]               # end params on the perimeter (α fixed)
        try
            out = findparscrit(S, xstar, p0, p1)
            pcrit = out.pars_crit      # (r1*, r2*, α*) or (Δr1*, Δr2*, α*) depending on setup
            # store Δr* only: convert r*→Δr* by subtracting r0
            push!(points, pcrit[1])# - r0[1])
            push!(points, pcrit[2])# - r0[2])
            push!(flags, out.flag)
        catch err
            @warn "Skipping ray due to error" δ1 δ2 err
        end
    end

    if !isempty(points)
        pts = reshape(points, 2, length(points)÷2)'

        # Masks per category (from decide_boundary_type)
        m_boundary = flags .== :boundary
        m_negative = flags .== :negative
        m_nonconv  = flags .== :nonconverged
        m_complex  = flags .== :complex
        m_unknown  = flags .== :unknown

        if any(m_boundary)
            Plots.scatter!(p, pts[m_boundary,1], pts[m_boundary,2];
                           ms=4, mc=:black, ma=0.9, marker=:circle,
                           label="numeric: boundary")
        end
        if any(m_negative)
            Plots.scatter!(p, pts[m_negative,1], pts[m_negative,2];
                           ms=4, mc=:indianred, ma=0.9, marker=:utriangle,
                           label="numeric: negative")
        end
        if any(m_complex)
            Plots.scatter!(p, pts[m_complex,1], pts[m_complex,2];
                           ms=4, mc=:steelblue, ma=0.9, marker=:diamond,
                           label="numeric: complex")
        end
        if any(m_nonconv)
            Plots.scatter!(p, pts[m_nonconv,1], pts[m_nonconv,2];
                           ms=4, mc=:orange, ma=0.9, marker=:x,
                           label="numeric: nonconverged")
        end
        if any(m_unknown)
            Plots.scatter!(p, pts[m_unknown,1], pts[m_unknown,2];
                           ms=4, mc=:gray, ma=0.7, marker=:star5,
                           label="numeric: unknown")
        end
    end

    mkpath(dirname(outfile))
    Plots.savefig(p, outfile)
    @info "Saved perimeter boundary overlay" outfile
    return nothing
end

##############################
# NEW: quick setup/demo
##############################
# Set variables to sample parameters and create systems
# Sampling setup (2×2 case for boundary polynomials)
ℓ, m, n = 2, 2, 2
seed    = 1
lower, upper = -1.0, 1.0

# Sample one parameter set
r0, A, B = generate_parameter_set(ℓ, m, n; lower=lower, upper=upper, seed=seed)
r0 = [-2., -1.0]   # override r0 for nicer plotting
#A (2×2)
A = [ 2        4/5
      1/2      1/8 ]

# B (2×2×2) — unspecified entries are 0
B = zeros(2,2,2)
B[1,1,1] = 4/3
B[1,2,1] = 1
B[2,2,1] = 10/7
B[1,1,2] = 1/3
B[1,2,2] = 4/3
B[2,2,2] = 9/4

save_boundary_overlay_perimeter(r0, A, B;
    α=0.1, halfwidth=5.0, n_per_side=100,
    outfile="plots/boundary_perimeter_points.pdf")


# # Build parameterized ODE polynomials f(x; Δr, α)
# @var x[1:n] Δr[1:n] a
# num_eqs = build_num_eqs(r0, A, B, x)
# S = get_parametrized_system(num_eqs, x, Δr, a)

# # Compute equilibria at Δr = 0 and chosen α
# α0 = 0.1
# sols, _ = compute_equilibria(r0, A, B; α=α0, real_only=true)
# @info "Real equilibria at α=$α0" sols

# # # Plot boundary curves in (Δr₁, Δr₂) analytically
# plot_boundaries_implicit(r0, A, B; α=α0, xlims=(-2,2), ylims=(-2,2))

#Calculate boundary numerically
#sample perturbations on a box centered around 0 of size 4 (xlims=(-2,2), ylims=(-2,2))
#calculate rcrit for each perturbation
#plot boundary points on top of implicit plot
#use similar code to this.
        # for idx in ind_bound
        #     Δr     = _embed2D_to_n(perts2D[idx], n)
        #     target = initial_parameters .+ Δr
        #     out    = U.findparscrit(syst, start_sol, initial_parameters, target)
        #     push!(parscrit_points, out.pars_crit)
        #     push!(flags, string(out.flag))
        # end

# save_boundaries_gif_implicit(r0, A, B; alphas=0.1:0.1:0.9,
#                              xlims=(-2,2), ylims=(-2,2),
#                              outfile="plots/boundaries.gif", fps=2)

#                              using ImplicitPlots, Plots, Random

# # --- batch runner ---
# """
#     batch_save_boundary_gifs(; seeds=1:100,
#                               α_range = range(0.01, 0.99; length=100),
#                               fps = 4,
#                               xlims = (-2.0, 2.0),
#                               ylims = (-2.0, 2.0),
#                               outdir = "plots",
#                               lower = -1.0, upper = 1.0)

# For each seed in `seeds`:
#   1) sample (r0, A, B) with generate_parameter_set(ℓ=2,m=2,n=2),
#   2) render a GIF of both boundaries across α ∈ α_range at `fps`,
#   3) save to "<outdir>/seed_<###>/boundaries.gif".
# """
# function batch_save_boundary_gifs(; seeds=1:100,
#                                    α_range = range(0.01, 0.99; length=100),
#                                    fps::Int = 4,
#                                    xlims = (-2.0, 2.0),
#                                    ylims = (-2.0, 2.0),
#                                    outdir::AbstractString = "boundary_gifs",
#                                    lower::Float64 = -1.0,
#                                    upper::Float64 = 1.0)
#     default(size=(720, 540), legend=:topright, framestyle=:box)
#     mkpath(outdir)

#     for s in seeds
#         try
#             println("→ seed $s")
#             # sample params (2×2 to match J1/J2)
#             ℓ = 2; m = 2; n = 2
#             r0, A, B = generate_parameter_set(ℓ, m, n; lower=lower, upper=upper, seed=s)

#             # per-seed dir + output gif
#             outfile = "$(outdir)/seed_$(s)_boundaries.gif"
#             save_boundaries_gif_implicit(r0, A, B;
#                 alphas  = α_range,
#                 xlims   = xlims,
#                 ylims   = ylims,
#                 outfile = outfile,
#                 fps     = fps
#             )
#             println("✓ saved $outfile")
#         catch err
#             @warn "Skipping seed due to error" seed=s error=err
#             # optionally show stacktrace:
#             # @warn sprint(showerror, err, catch_backtrace())
#             continue
#         end
#     end
#     return nothing
# end


# # --- run it ---
# batch_save_boundary_gifs(; seeds=1:100,
#                           α_range=range(0.01, 0.99; length=100),
#                           fps=10,
#                           xlims=(-2, 2),
#                           ylims=(-2, 2),
#                           outdir="boundary_gifs")
