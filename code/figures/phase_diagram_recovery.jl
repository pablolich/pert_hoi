using HomotopyContinuation
using Combinatorics, Plots

# === STEP 1: Build square mesh ===
function build_square_mesh(npoints::Float64, boxsize::Float64)
    lengthvec = Int(round(sqrt(npoints)))
    outer_vec = collect(range(-boxsize / 2, boxsize / 2, length = lengthvec))
    perts = Array{Tuple{Float64, Float64}}(undef, lengthvec, lengthvec)
    for i in 1:lengthvec, j in 1:lengthvec
        perts[i, j] = (outer_vec[i], outer_vec[j])
    end
    return perts, outer_vec
end

# === STEP 2: Greedy reordering based on neighbors ===
function recover_phase_diagram(V::AbstractMatrix{<:AbstractVector{<:Real}})
    N = size(V, 1)
    nsol = length(V[1,1])
    perms = collect(Combinatorics.permutations(1:nsol))
    M = [zeros(Float64, N, N) for _ in 1:nsol]

    for k in 1:nsol
        M[k][1,1] = V[1,1][k]
    end

    for i in 1:N, j in 1:N
        if i == 1 && j == 1
            continue
        end
        reference = if i > 1
            [M[k][i-1,j] for k in 1:nsol]
        else
            [M[k][i,j-1] for k in 1:nsol]
        end

        values = V[i,j]
        best_perm, best_cost = nothing, Inf
        for perm in perms
            permuted = values[perm]
            cost = sum(abs.(permuted .- reference))
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

# === STEP 3: Setup symbolic system ===
n = 2
d = 2
@var x[1:n] α[1:d]
ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n, coeff_name = :c)
rng = MersenneTwister(2)
pars = sample_parameters(n, d+1, rng)
eqs = build_glvhoi(pars, x)

inds_growth = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_param = parametrize_interactions(eqs, ref_eqs, x, inds_growth)
eqs_param_str = parametrize_strengths(eqs_param, x, α)
syst = build_parametrized_glvhoi(eqs_param_str, x, coeffs_mat, inds_growth, α)

alpha = 0.2
syst = evaluate_pars(syst, α, [1 - alpha, alpha])
initial_parameters = pars[1]

# === STEP 4: Solve across perturbation grid and extract Re(x₁) ===
nperts = 40.0^2
pert_size = 3.0
perts, outer_vec = build_square_mesh(nperts, pert_size)
lengthvec = size(perts, 1)

#Now obtain the boundary using findparscrit

#get boundary of perturbations
ind_bound = [i for i in 1:Int(nperts) if any(abs.(perts[i]) .== pert_size/2)]
bound_perts = perts[ind_bound]
start_solution = [1,1]

parscrit_points = []
flags = []

for i in 1:length(ind_bound)
    println("Boundary perturbation: ", i)
    end_parameters = pars[1] .+ bound_perts[i]
    pars_crit, xstar_crit, flag = findparscrit(syst, start_solution, initial_parameters, end_parameters)
    push!(parscrit_points, pars_crit)
    push!(flags, string(flag))  # convert Symbol to String for plotting
end

using DataStructures

grouped_points = DefaultDict{String, Vector{Tuple{Float64, Float64}}}(() -> [])

for (pt, f) in zip(parscrit_points, flags)
    push!(grouped_points[f], (pt[1], pt[2]))
end

# Prepare 4 fields: real/imag of species 1 and 2 for solution 1
V_real = [Array{Float64}(undef, lengthvec, lengthvec) for _ in 1:2]  # species 1, 2
V_imag = [Array{Float64}(undef, lengthvec, lengthvec) for _ in 1:2]

for i in 1:lengthvec, j in 1:lengthvec
    println("Solving at ($i, $j)...")
    r1, r2 = perts[i, j]
    perturbed_pars = initial_parameters .+ [r1, r2]
    syst_fixed = fix_parameters(syst, perturbed_pars)
    sols = [s.solution for s in path_results(solve(syst_fixed))]

    sol1 = sols[1]  # only solution 1
    for s in 1:2
        V_real[s][i, j] = real(sol1[s])
        V_imag[s][i, j] = imag(sol1[s])
    end
end


# === STEP 5: Recover 4 smooth phase diagrams ===
M1_real_fields = recover_phase_diagram(V_real_species1)

# === STEP 6: Plot or save each recovered field ===
for p in 1:4
    p_h = heatmap(
        outer_vec .+ initial_parameters[1],
        outer_vec .+ initial_parameters[2],
        M1_real_fields[p]',
        title = "Re(x₁), solution $p",
        xlabel = "r₁",
        ylabel = "r₂",
        color = cgrad(:viridis),
        clims = extrema(M1_real_fields[p])
    )
    for (label, pts) in grouped_points
        scatter!(
            p_h,
            first.(pts),
            last.(pts),
            label = label,
            markersize = 4,
            alpha = 0.8,
        )
    end
    scatter!(
        p_h, [pars[1][1]], [pars[1][2]], 
        label = "Origin", 
        color = :black, markershape = :x, markersize = 6)
    
    savefig("Re_x1_solution_$p.png")
end
