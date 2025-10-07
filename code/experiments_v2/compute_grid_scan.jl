# compute_grid_scan.jl
#
# For each α in `alphavec`, over a square grid of perturbations:
#   - Track the chosen unperturbed equilibrium via parameter homotopy
#   - Record min(real.(x)) and min(abs.(imag.(x))) at the target parameters
#   - On the grid boundary, compute feasibility boundary with U.findparscrit
# Saves a single serialized artifact per (n, seed) containing all α slices.

include("src_functions.jl")

using HomotopyContinuation
using HomotopyContinuation: path_results
using Serialization, LinearAlgebra, Statistics, Random

# -------------------------
# user-configurable settings
# -------------------------
const ORDER_D      = 2                    # polynomial order (degree)
const D_PERT       = 0                    # perturb growth rates
const N_SIDE       = 40                   # grid side length => N_SIDE^2 total perturbations
const PERT_SIZE    = 10.                   # box size edge length
const ALPHAVEC     = collect(range(0, 0.3; length=10))  # α sweep
const SEED_ID      = 3                    # which (r,A,B,…) seed file to use
const OUTDIR       = "grid_scans"

# -------------------------
# helpers
# -------------------------

# Choose a start solution at the unperturbed parameters: the one closest to 1⃗
# (works well with your construction where x*≈1⃗ is planted).
function _pick_start_solution(sols::Vector{<:AbstractVector{<:Complex}}, n::Int)
    dists = [norm(real.(s) .- ones(n)) for s in sols]
    idx   = argmin(dists)
    return sols[idx]
end

# Build an n-dimensional perturbation from a 2D grid point (only first 2 coords vary)
@inline function _embed2D_to_n(δ::NTuple{2,Float64}, n::Int)
    p = zeros(Float64, n)
    p[1] = δ[1]; p[2] = δ[2]
    return p
end

# Boundary test for 2D grid point on the square (half-box at ±PERT_SIZE/2)
@inline _on_box_boundary(δ::NTuple{2,Float64}, half::Float64) =
    isapprox(abs(δ[1]), half; atol=1e-9) || isapprox(abs(δ[2]), half; atol=1e-9)

# -------------------------
# main driver
# -------------------------
function compute_grid_scan(n::Int; seed::Int=SEED_ID,
                           n_side::Int=N_SIDE, pert_size::Float64=PERT_SIZE,
                           alphavec::Vector{Float64}=ALPHAVEC,
                           order_d::Int=ORDER_D, d_pert::Int=D_PERT)

    # 0) Load (or ensure) parameter tuple (r,A,B,...) for this n
    param_dir = joinpath(@__DIR__, "parameter_sets", "n_$(n)")
    isdir(param_dir) || mkpath(param_dir)
    # if you want to regenerate, uncomment:
    # U.generate_and_save_parameter_sets(seed, n, n, n, param_dir; lower=-1.0, upper=1.0, prefix="seed")

    param_path = joinpath(param_dir, "seed_seed$(seed).bin")
    (r, A, B) = open(param_path, "r") do io
        deserialize(io)
    end

    # 1) Build symbolic system data once for this n
    @var x[1:n]
    @var α[1:2]
    ref_eqs, coeffs_mat = U.get_ref_polynomials(x, order_d, n; coeff_name=:c)
    num_eqs             = U.build_glvhoi((r, A, B), x)
    syst_generic        = U.get_parametrized_system(num_eqs, ref_eqs, coeffs_mat,
                                                    d_pert, x, α)

    # 2) Prepare the 2D grid of perturbations (centered at 0)
    # You said to assume U.build_square_grid exists; it should return a Vector{NTuple{2,Float64}}
    npoints  = n_side^2
    perts2D  = U.build_square_mesh(npoints, pert_size)  # assumed available in src_functions.jl
    halfbox  = pert_size / 2
    outervec = collect(range(-halfbox, halfbox; length=n_side)) # for axes

    # Index of boundary points on the square
    ind_bound = findall(i -> _on_box_boundary(perts2D[i], halfbox), eachindex(perts2D))

    # 3) Sweep α; collect slices
    zreal_slices = Vector{Matrix{Float64}}()
    zimag_slices = Vector{Matrix{Float64}}()
    boundary_pts_slices   = Vector{Matrix{Float64}}()  # each is (#boundary)×2 (first two dims)
    boundary_flags_slices = Vector{Vector{String}}()

    for αval in alphavec
        println("Computing α = $αval for n = $n, seed = $seed...")
        # Fix α strengths in the system
        syst = U.evaluate_pars(syst_generic, α, [1-αval, αval])

        # Unperturbed starting solutions at parameters=r
        S_fixed   = U.plug_parameters_into_system(syst, r)
        res0      = solve(S_fixed)
        sols0     = solutions(res0)
        start_sol = _pick_start_solution(sols0, n)
        initial_parameters = r

        # Track over every grid perturbation
        min_real  = Vector{Float64}(undef, npoints)
        min_imag  = Vector{Float64}(undef, npoints)

        for i in 1:npoints
            Δr      = _embed2D_to_n(perts2D[i], n)
            target  = initial_parameters .+ Δr
            # Follow the chosen unperturbed equilibrium via coefficient homotopy
            respaths = path_results(solve(syst, [start_sol];
                                          start_parameters = initial_parameters,
                                          target_parameters = target))
            # first (and only) path result corresponds to start_sol
            sol = respaths[1].solution
            # a) min component of the real part
            min_real[i] = minimum(real.(sol))
            # b) (requested) "minimum imaginary part": take min(abs(imag)) across components
            min_imag[i] = minimum(abs.(imag.(sol)))
        end

        # b) boundary via findparscrit only along the square boundary rays
        parscrit_points = Vector{Vector{Float64}}()
        flags           = String[]
        for idx in ind_bound
            Δr     = _embed2D_to_n(perts2D[idx], n)
            target = initial_parameters .+ Δr
            out    = U.findparscrit(syst, start_sol, initial_parameters, target)
            push!(parscrit_points, out.pars_crit)
            push!(flags, string(out.flag))
        end

        # reshape vectors into (y,x) matrices for heatmaps
        zreal = transpose(reshape(min_real, (n_side, n_side)))
        zimag = transpose(reshape(min_imag, (n_side, n_side)))

        push!(zreal_slices, zreal)
        push!(zimag_slices, zimag)

        # store only first 2 components for 2D overlay
        if !isempty(parscrit_points)
            pts = hcat(getindex.(parscrit_points, 1), getindex.(parscrit_points, 2))
            push!(boundary_pts_slices, pts)
        else
            push!(boundary_pts_slices, zeros(0, 2))
        end
        push!(boundary_flags_slices, flags)
    end

    # 4) Serialize artifact
    isdir(OUTDIR) || mkpath(OUTDIR)
    outpath = joinpath(OUTDIR, "grid_scan_n_$(n)_seed_$(seed).bin")
    open(outpath, "w") do io
        serialize(io, (
            n            = n,
            seed         = seed,
            order_d      = order_d,
            n_side       = n_side,
            pert_size    = pert_size,
            alphas       = alphavec,
            x_axis       = r[1] .+ outervec,
            y_axis       = r[2] .+ outervec,
            zreal_slices = zreal_slices,
            zimag_slices = zimag_slices,
            boundary_pts = boundary_pts_slices,
            boundary_flags = boundary_flags_slices,
        ))
    end
    println("✔ saved → ", outpath)
end

# ---------
# Example(s)
# ---------
compute_grid_scan(2)
