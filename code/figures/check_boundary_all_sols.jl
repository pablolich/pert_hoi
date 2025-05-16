using Plots, LinearAlgebra, DataStructures
include("../source_function/general_perturbation_functions.jl")

function build_square_mesh(npoints::Float64, boxsize::Float64)
    lengthvec = Int(round(sqrt(npoints)))
    outer_vec = collect(range(-boxsize/2, boxsize/2, length = lengthvec))
    return [collect(Iterators.product(outer_vec, outer_vec))...]
end

# === Setup system ===
n = 2
@var x[1:n]
d = 2
@var α[1:d]
ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n, coeff_name = :c)

rng = MersenneTwister(2)
pars = sample_parameters(n, d+1, rng)
eqs = build_glvhoi(pars, x)

inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
eqs_inter_str = parametrize_strengths(eqs_inter, x, α)
syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

alpha = 0.2
syst = evaluate_pars(syst, α, [1-alpha, alpha])
initial_parameters = pars[1]
start_solution = [1.0 + 0im, 1.0 + 0im]

# === Critical Boundary Points ===
nperts = 5.0^2
pert_size = 0.2
lengthvec = Int(round(sqrt(nperts)))
outer_vec = collect(range(-pert_size/2, pert_size/2, length = lengthvec))
perts = build_square_mesh(nperts, pert_size)
ind_bound = [i for i in 1:Int(nperts) if any(abs.(perts[i]) .== pert_size/2)]
bound_perts = perts[ind_bound]

parscrit_points = []
flags = []
for i in 1:length(ind_bound)
    println("Boundary perturbation: ", i)
    end_parameters = initial_parameters .+ bound_perts[i]
    pars_crit, xstar_crit, flag = findparscrit(syst, start_solution, initial_parameters, end_parameters)
    push!(parscrit_points, pars_crit)
    push!(flags, string(flag))
end

# === Heatmap via Aligned Solution Tracking ===

# Sort perturbations by Euclidean distance to initial_parameters
pert_distances = [norm(p) for p in perts]
sorted_inds = sortperm(pert_distances)
sorted_perts = perts[sorted_inds]

# Solve at unperturbed parameters and store all solutions
syst_fixed = fix_parameters(syst, initial_parameters)
initial_solutions = path_results(solve(syst_fixed))
previous = [s.solution for s in initial_solutions]

# Preallocate aligned solutions matrix (nperts × ns × nvars)
num_sols = length(previous)
aligned_solutions = Vector{Matrix{ComplexF64}}()

for (i, p) in enumerate(sorted_perts)
    println("Aligned solve: ", i)
    perturbed_pars = initial_parameters .+ p
    syst_fixed = fix_parameters(syst, perturbed_pars)
    current_solutions = [s.solution for s in path_results(solve(syst_fixed))]
    println("Previous solutions: ", previous)
    println("Current solutions: ", current_solutions)
    println("Perturbation: ", p)

    # Align current to previous using minimum pairwise distances
    D = [norm(current - prev) for current in current_solutions, prev in previous]
    assignment = argmin.(eachrow(D))  # which previous each current is closest to

    # Now reorder `current_solutions` to align to `previous`
    aligned = similar(previous)
    for j in 1:num_sols
        # Find best match for previous[j]
        idx = argmin([norm(current_solutions[k] - previous[j]) for k in 1:num_sols])
        aligned[j] = current_solutions[idx]
    end
    # `aligned` is a Vector of 4 vectors (length 2 each)
    # Turn into a 2 × 4 matrix
    aligned_matrix = hcat(aligned...)  # shape: (2, 4)

    # Store it
    push!(aligned_solutions, aligned_matrix)
end

# `aligned_solutions_sorted` was built in sorted_perts order
# Reorder them back to original perts grid layout:
aligned_solutions_grid = similar(aligned_solutions)
for (i_sorted, i_grid) in enumerate(sorted_inds)
    aligned_solutions_grid[i_grid] = aligned_solutions[i_sorted]
end
