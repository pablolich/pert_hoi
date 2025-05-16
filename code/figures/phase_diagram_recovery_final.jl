# ========== SETUP ========== #
using HomotopyContinuation, Combinatorics
using DataStructures
using CairoMakie

const cm = CairoMakie

# Build square mesh
function build_square_mesh(npoints::Float64, boxsize::Float64)
    lengthvec = Int(round(sqrt(npoints)))
    outer_vec = collect(range(-boxsize / 2, boxsize / 2, length = lengthvec))
    perts = Array{Tuple{Float64, Float64}}(undef, lengthvec, lengthvec)
    for i in 1:lengthvec, j in 1:lengthvec
        perts[i, j] = (outer_vec[i], outer_vec[j])
    end
    return perts, outer_vec
end

using Combinatorics

using Combinatorics

function recover_phase_diagram_twopass(V::AbstractMatrix{<:AbstractVector{<:Real}})
    N = size(V, 1)
    nsol = length(V[1,1])
    perms = collect(Combinatorics.permutations(1:nsol))

    # Helper to compute cost against available neighbors
    function compute_cost(i, j, permuted, M)
        neighbors = []
        if i > 1
            push!(neighbors, [M[k][i-1,j] for k in 1:nsol])
        end
        if j > 1
            push!(neighbors, [M[k][i,j-1] for k in 1:nsol])
        end
        if i < N
            push!(neighbors, [M[k][i+1,j] for k in 1:nsol])
        end
        if j < N
            push!(neighbors, [M[k][i,j+1] for k in 1:nsol])
        end
        # Cost = sum of L2 distances to neighbors
        return sum([sum(abs2, permuted .- neighbor) for neighbor in neighbors])
    end

    function pass_align(V, direction=:forward)
        M = [zeros(Float64, N, N) for _ in 1:nsol]
        visited = falses(N, N)

        indices = direction == :forward ? ((i,j) for i in 1:N, j in 1:N) :
                                          ((i,j) for i in N:-1:1, j in N:-1:1)

        for (i, j) in indices
            values = V[i,j]
            if i == 1 && j == 1 && direction == :forward
                for k in 1:nsol
                    M[k][i,j] = values[k]
                end
                continue
            elseif i == N && j == N && direction == :reverse
                for k in 1:nsol
                    M[k][i,j] = values[k]
                end
                continue
            end

            best_perm, best_cost = nothing, Inf
            for perm in perms
                permuted = values[perm]
                cost = compute_cost(i, j, permuted, M)
                if cost < best_cost
                    best_cost = cost
                    best_perm = perm
                end
            end
            for k in 1:nsol
                M[k][i,j] = values[best_perm[k]]
            end
        end
        return M
    end

    # Run both passes
    M_fwd = pass_align(V, :forward)
    M_rev = pass_align(V, :reverse)

    # Final decision: pick lower-cost result at each point
    M_final = [zeros(Float64, N, N) for _ in 1:nsol]
    for i in 1:N, j in 1:N
        v = V[i,j]
        v_fwd = [M_fwd[k][i,j] for k in 1:nsol]
        v_rev = [M_rev[k][i,j] for k in 1:nsol]

        cost_fwd = compute_cost(i, j, v_fwd, M_fwd)
        cost_rev = compute_cost(i, j, v_rev, M_rev)

        chosen = cost_fwd ≤ cost_rev ? v_fwd : v_rev
        for k in 1:nsol
            M_final[k][i,j] = chosen[k]
        end
    end

    return M_final
end

function matrices_to_vectorgrid(M::Vector{Matrix{Float64}})
    N = size(M[1], 1)
    V = Matrix{Vector{Float64}}(undef, N, N)
    for i in 1:N, j in 1:N
        V[i,j] = [M[k][i,j] for k in 1:length(M)]
    end
    return V
end

function multipass_alignment(V_init::Vector{Matrix{Vector{Float64}}}, npasses::Int = 4)
    V = V_init
    M = nothing  # initialize M outside the loop

    for pass in 1:npasses
        M = [recover_phase_diagram_twopass(V[s]) for s in 1:2]
        V = [matrices_to_vectorgrid(M[s]) for s in 1:2]
    end

    return M
end




# function recover_phase_diagram(V::AbstractMatrix{<:AbstractVector{<:Real}})
#     N = size(V, 1)
#     nsol = length(V[1,1])
#     perms = collect(Combinatorics.permutations(1:nsol))
#     M = [zeros(Float64, N, N) for _ in 1:nsol]

#     for k in 1:nsol
#         M[k][1,1] = V[1,1][k]
#     end

#     for i in 1:N, j in 1:N
#         if i == 1 && j == 1
#             continue
#         end

#         reference = if i > 1 && j > 1
#             0.5 .* ([M[k][i-1,j] for k in 1:nsol] .+ [M[k][i,j-1] for k in 1:nsol])
#         elseif i > 1
#             [M[k][i-1,j] for k in 1:nsol]
#         elseif j > 1
#             [M[k][i,j-1] for k in 1:nsol]
#         end

#         values = V[i,j]
#         best_perm, best_cost = nothing, Inf
#         for perm in perms
#             permuted = values[perm]
#             cost = sum(abs2, permuted .- reference)
#             if cost < best_cost
#                 best_cost = cost
#                 best_perm = perm
#             end
#         end

#         for k in 1:nsol
#             M[k][i,j] = values[best_perm[k]]
#         end
#     end

#     return M
# end

# ========== SYSTEM DEFINITION ========== #
n, d = 2, 2
@var x[1:n] 
@var α[1:d]
ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n, coeff_name = :c)
rng = MersenneTwister(5)
pars = sample_parameters(n, d+1, rng)
eqs = build_glvhoi(pars, x)
inds_growth = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_param = parametrize_interactions(eqs, ref_eqs, x, inds_growth)
eqs_param_str = parametrize_strengths(eqs_param, x, α)
syst = build_parametrized_glvhoi(eqs_param_str, x, coeffs_mat, inds_growth, α)

alpha = 0.2
syst = evaluate_pars(syst, α, [1 - alpha, alpha])
initial_parameters = pars[1]

# ========== GRID SETUP ========== #
nperts = 40.0^2
pert_size = 4.0
perts, outer_vec = build_square_mesh(nperts, pert_size)
lengthvec = size(perts, 1)

# --- Compute critical boundaries ---
ind_bound = [i for i in 1:lengthvec^2 if any(abs.(Tuple(perts[i])) .== pert_size/2)]
bound_perts = [perts[i] for i in ind_bound]
start_solution = [1.0, 1.0]

parscrit_points = []
flags = []

for i in 1:length(bound_perts)
    println("Boundary perturbation: ", i)
    end_parameters = initial_parameters .+ collect(bound_perts[i])
    pars_crit, _, flag = findparscrit(syst, start_solution, initial_parameters, end_parameters)
    push!(parscrit_points, pars_crit)
    push!(flags, string(flag))
end

grouped_points = DefaultDict{String, Vector{Tuple{Float64, Float64}}}(() -> [])
for (pt, f) in zip(parscrit_points, flags)
    push!(grouped_points[f], (pt[1], pt[2]))
end

# ========== SOLVE & EXTRACT ========== #
V_real = [Array{Vector{Float64}}(undef, lengthvec, lengthvec) for _ in 1:2]
V_imag = [Array{Vector{Float64}}(undef, lengthvec, lengthvec) for _ in 1:2]

for i in 1:lengthvec, j in 1:lengthvec
    r1, r2 = perts[i, j]
    perturbed_pars = initial_parameters .+ [r1, r2]
    syst_fixed = fix_parameters(syst, perturbed_pars)
    sols = [s.solution for s in path_results(solve(syst_fixed))]

    for s in 1:2
        V_real[s][i, j] = [real(sol[s]) for sol in sols]
        V_imag[s][i, j] = [imag(sol[s]) for sol in sols]
    end
end

# ========== RECOVER FIELDS & PLOT ========== #
x_axis = outer_vec .+ initial_parameters[1]
y_axis = outer_vec .+ initial_parameters[2]
titles = ["Re(x1)", "Im(x1)", "Re(x2)", "Im(x2)"]

M_real = multipass_alignment(V_real, 10)
M_imag = multipass_alignment(V_imag, 10)

fig = CairoMakie.Figure(resolution = (4000, 1000))
main_layout = fig[1, 1] = CairoMakie.GridLayout()

titles = ["Re(x1)", "Im(x1)", "Re(x2)", "Im(x2)"]

for sol_index in 1:4
    col = sol_index
    outer_grid = main_layout[1, col] = CairoMakie.GridLayout()
    grid = outer_grid[1, 1] = CairoMakie.GridLayout()

    hms = Heatmap[]
    for (i, Z) in enumerate([M_real[1][sol_index], M_imag[1][sol_index], M_real[2][sol_index], M_imag[2][sol_index]])
        ax = CairoMakie.Axis(grid[i <= 2 ? 1 : 2, i % 2 == 1 ? 1 : 2],
                  title = titles[i],
                  xlabel = "r₁", ylabel = "r₂")
        hm = CairoMakie.heatmap!(ax, x_axis, y_axis, Z, colormap = :inferno)
        push!(hms, hm)
        for (label, pts) in grouped_points
            CairoMakie.scatter!(ax, first.(pts), last.(pts), color = label == "boundary" ? :black : :red, markersize = 5)
        end
        CairoMakie.scatter!(ax, [initial_parameters[1]], [initial_parameters[2]], color = :black, marker = :x, markersize = 10)
    end

    # Add a colorbar below the 2x2 grid (use first heatmap handle directly)
    CairoMakie.Colorbar(outer_grid[2, 1], hms[1], vertical = false, height = 30)
end

CairoMakie.save("all_solutions_2x2_grids.png", fig)
