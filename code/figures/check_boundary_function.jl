
#CHECK FUNCTION TO CALCULATE FEASIBILITY BOUNDARY AND PLOT IT
include("../source_function/general_perturbation_functions.jl")
using Plots 

"""
create a 2D square mesh around 0 of size boxsize containing sqrt(npoints)
"""
function build_square_mesh(npoints::Float64, boxsize::Float64)
    lengthvec = Int(round(sqrt(npoints)))
    outer_vec = collect(range(-boxsize/2,boxsize/2, length = lengthvec))
    return [collect(Iterators.product(outer_vec, outer_vec))...]
end

#for each point in the box, track solution from unperturbed to that point
#test all functions 

#sample model parameters
n = 2
@var  x[1:n]
d = 2
@var α[1:d]
#get reference coefficients for systems of such dimension & degree
ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n, coeff_name = :c)

#set a random seed
rng = MersenneTwister(5) #6, 7
#seed=2 with pert_size=1.0 generates only positive solutions
using Random

rng = MersenneTwister(1234)  # or any desired seed
ϵ = 0.

r = [-5., -5.]

A = [
    2.0    4/5;
    1/2    1/8
] .+ ϵ .* randn(rng, 2, 2)

B = zeros(2, 2, 2)
B[1,1,1] = 4/3
B[2,1,1] = 1.0
B[2,2,1] = 10/7
B[1,1,2] = 1/3
B[2,1,2] = 4/3
B[2,2,2] = 9/4
B .+= ϵ .* randn(rng, 2, 2, 2)

pars = (r, A, B)
pars_rand = sample_parameters(n, d+1, rng)


#build glv system with sampled parameters
eqs = build_glvhoi(pars, x)

#build parametrized system through combination of symbols and numbers
inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
eqs_inter_str = parametrize_strengths(eqs_inter, x, α)
syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

#sample square grid of perturbations
nperts = 40.0^2
pert_size = 10.5
lengthvec = Int(round(sqrt(nperts)))
outer_vec = collect(range(-pert_size/2, pert_size/2, length = lengthvec))
perts = build_square_mesh(nperts, pert_size)

#evaluate model at a certain HOI strength
alpha = 0.1
syst = evaluate_pars(syst, α, [1-alpha, alpha])
#initialize vector of solutions
min_comps_real  = []
min_comps_imag = []

#start a particular solution (we planted the 1,1 solution) and parameters
#solve to get initial solutions
initial_parameters = pars[1]
# Fix the parameters in the system
S_fixed = fix_parameters(syst, initial_parameters)

# Solve the system with fixed parameters
result = solve(S_fixed)
start_solutions = [[0.11192548903069875 + 0.0im, 4.401106963819901 + 3.851859888774472e-34im]]

for i in 1:Int(nperts)
    println("Perturbation: ", i)
    #specify end parameters
    end_parameters = pars[1] .+ perts[i]
    res = path_results(
        solve(
            syst, start_solutions;
            start_parameters = initial_parameters,  
            target_parameters = end_parameters
        )
    )[1]
    #select the minimum component
    push!(min_comps_real, minimum(real(res.solution)))
    push!(min_comps_imag, maximum(abs.(imag(res.solution))))
end

#Now obtain the boundary using findparscrit

#get boundary of perturbations
ind_bound = [i for i in 1:Int(nperts) if any(abs.(perts[i]) .== pert_size/2)]
bound_perts = perts[ind_bound]
start_solution = [0.11192548903069875 + 0.0im, 4.401106963819901 + 3.851859888774472e-34im]
#start_solution = [1,1]

parscrit_points = []
flags = []

for i in 1:length(ind_bound)
    println("Boundary perturbation: ", i)
    end_parameters = pars[1] .+ bound_perts[i]
    pars_crit, xstar_crit, flag = findparscrit(syst, start_solution, initial_parameters, end_parameters)
    push!(parscrit_points, pars_crit)
    push!(flags, string(flag))  # convert Symbol to String for plotting
end


lengthvec = Int(round(sqrt(nperts)))

# build the axes
x_axis = pars[1][1] .+ outer_vec
y_axis = pars[1][2] .+ outer_vec

# reshape min_comps into a square matrix
zreal = reshape(min_comps_real, (lengthvec, lengthvec))
zreal = transpose(zreal)  # heatmap expects (y, x) layout
zimag = reshape(min_comps_imag, (lengthvec, lengthvec))
zimag = transpose(zimag)  # heatmap expects (y, x) layout

# Determine symmetric color limits around 0
zminreal = minimum(zreal)
zmaxreal = maximum(zreal)
absmaxreal = max(abs(zminreal), abs(zmaxreal))

# Create the heatmap with symmetric colorbar centered at 0
p1 = heatmap(
    x_axis,
    y_axis,
    zreal,
    c = cgrad(:RdBu, rev = true),
    clims = (-absmaxreal, absmaxreal)  # center color at 0
)

# Determine symmetric color limits around 0
zminimag = minimum(zimag)
zmaximag = maximum(zimag)
absmaximag = max(abs(zminimag), abs(zmaximag))

flag_labels = string.(flags)

# Group points by flag
using DataStructures

# Create a default dictionary where each missing key maps to an empty vector of (x, y) tuples
grouped_points = DefaultDict{String, Vector{Tuple{Float64, Float64}}}(() -> Vector{Tuple{Float64, Float64}}())

# Fill the dictionary by grouping by flag (converted to string)
for (pt, flag) in zip(parscrit_points, flags)
    push!(grouped_points[string(flag)], (pt[1], pt[2]))
end

# Plot each group separately with a distinct legend label
for (label, pts) in grouped_points
    scatter!(
        p1,
        first.(pts),
        last.(pts),
        label = label,
        markersize = 4
    )
 
end

scatter!(p1, [pars[1][1]], [pars[1][2]], label = "Origin", color = :black, markershape = :x, markersize = 6)

plot(p1)
savefig("boundary_real_imag.png")
# Save the figure to a file
#savefig(p, "boundary_detection_plot.png")

# #plot
# heatmap(pars[1][1] .+ outer_vec, pars[1][2] .+ outer_vec, 
# transpose(reshape(min_comps, Int(sqrt(nperts)), Int(sqrt(nperts)))),
# c =cgrad(:RdBu, rev = true))
# scatter!(parscrit_points[:,1], parscrit_points[:,2], c=:black)
