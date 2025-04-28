
#CHECK FUNCTION TO CALCULATE FEASIBILITY BOUNDARY AND PLOT IT
include("general_perturbation_functions.jl")
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
rng = MersenneTwister(1)
#seed=2 with pert_size=1.0 generates only positive solutions
pars = sample_parameters(n, d+1, rng)
#build glv system with sampled parameters
eqs = build_glvhoi(pars, x)

#build parametrized system through combination of symbols and numbers
inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
eqs_inter_str = parametrize_strengths(eqs_inter, x, α)
syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

#sample square grid of perturbations
nperts = 40.0^2
pert_size = 10.0
lengthvec = Int(round(sqrt(nperts)))
outer_vec = collect(range(-pert_size/2, pert_size/2, length = lengthvec))
perts = build_square_mesh(nperts, pert_size)

#evaluate model at a certain HOI strength
alpha = 0.5
syst = evaluate_pars(syst, α, [1-alpha, alpha])
#initialize vector of solutions
min_comps = []

#start a particular solution (we planted the 1,1 solution) and parameters
start_solutions = [[1.0 + 0.0im, 1.0 + 0.0im]]
initial_parameters = pars[1]

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
    push!(min_comps, minimum(real(res.solution)))
end

#Now obtain the boundary using findparscrit

#get boundary of perturbations
ind_bound = [i for i in 1:Int(nperts) if any(abs.(perts[i]) .== pert_size/2)]
bound_perts = perts[ind_bound]
start_solution = [1,1]

parscrit_points = Vector{Float64}[]
for i in 1:length(ind_bound)
    println("Boundary perturbation: ", i)
    end_parameters = pars[1] .+ bound_perts[i]
    pars_crit, xstar_crit, flag = findparscrit(syst, start_solution, initial_parameters, end_parameters)
    push!(parscrit_points, pars_crit)
end

lengthvec = Int(round(sqrt(nperts)))

# build the axes
x_axis = pars[1][1] .+ outer_vec
y_axis = pars[1][2] .+ outer_vec

# reshape min_comps into a square matrix
z = reshape(min_comps, (lengthvec, lengthvec))
z = transpose(z)  # heatmap expects (y, x) layout

# Create the heatmap
p = heatmap(
    x_axis,
    y_axis,
    z,
    c = cgrad(:RdBu, rev=true)
)

# Overlay critical points
scatter!(
    p,
    first.(parscrit_points),   # x coordinates
    last.(parscrit_points),    # y coordinates
    c = :black,
    markersize = 3
)

# Save the figure to a file
savefig(p, "boundary_detection_plot.png")

# #plot
# heatmap(pars[1][1] .+ outer_vec, pars[1][2] .+ outer_vec, 
# transpose(reshape(min_comps, Int(sqrt(nperts)), Int(sqrt(nperts)))),
# c =cgrad(:RdBu, rev = true))
# scatter!(parscrit_points[:,1], parscrit_points[:,2], c=:black)
