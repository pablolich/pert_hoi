#this script deals with the question: 
#how does the histogram of distances change when
#d = 2, we perturb growth rates (d_pert = 0), set n = 2, the maximum perturbation magnitude is 10
#and vary alpha, perturbation direction, and look at different parameter sets, indexed by their seed.

using DelimitedFiles

include("../source_function/general_perturbation_functions.jl")

"""
wrapper for readibility; runs in series all the functions necessary to prepare the system for search
"""
function get_parametrized_system(
    num_eqs::Vector{Expression}, 
    symb_eqs::Vector{Expression}, 
    coeffs_mat,
    perturb_order::Int64, 
    vars::Vector{Variable},
    strengths::Vector{Variable})
    #get indices of growth rates
    inds_growth_rates = get_ind_coeffs_subs(symb_eqs[1], vars, "order", [perturb_order])
    #parametrize growth rates so we can perturb them
    eqs_inter = parametrize_interactions(num_eqs, symb_eqs, vars, inds_growth_rates)
    #parametrize relative strength of interaction orders so we can vary them
    eqs_inter_str = parametrize_strengths(eqs_inter, vars, strengths)
    #return System of equations
    return build_parametrized_glvhoi(eqs_inter_str, vars, coeffs_mat, inds_growth_rates, strengths)
end

#set fix parameters for simulations
nsim = 1000 #number of simulations
d = 2; @var α[1:d] #degree of polynomials to solve; corresponding strengths
n = 5; @var x[1:n] #number of equations; corresponding variables
pert_size = 10.0 #maximum perturbation
d_pert = 0 #order of parameters to perturb
n_perts = 1000 #number of perturbations
pert_dirs = points_hypersphere(n, 1.0, n_perts, false) #perturbation directions
alpha_vec = [0.1, 0.9] #values of relative interaction strength
init_sol = repeat([1], n) #set initial solution to 1

#build skeleton polynomial system for n, d
ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n)

open("../../data/results_simplest_extension_n_5.csv", "a") do io
    #run all simulations
    for seed_i in 1:nsim
        rng = MersenneTwister(seed_i)
        #sample parameters. note that they are deterministic given n, order and rng.
        pars = sample_parameters(n, d+1, rng)
        initial_pars = pars[1]
        #build numerical glvhoi system
        eqs = build_glvhoi(pars, x)

        #build system parametrizing strengths and parameters to be perturbed
        syst = get_parametrized_system(eqs, ref_eqs, coeffs_mat, d_pert, x, α)

        for pert_i in 1:n_perts 
            end_parameters = initial_pars .+ pert_size .* pert_dirs[pert_i,:]

            for alpha_i in alpha_vec #loop through each relative strength value
                syst_alpha = evaluate_pars(syst, α, [1-alpha_i, alpha_i])
                println("Simulation: ", seed_i, " parameter set: ", pert_i, " relative strength: ", alpha_i)
                pars_crit = findparscrit(syst_alpha, init_sol, initial_pars, end_parameters)
                if pars_crit == -1
                    return pars, pert_size, pert_dirs[pert_i,:], alpha_i
                end
                δₚ = norm(pars_crit .- initial_pars)

                iteration_result = [seed_i n alpha_i pert_i δₚ]
                writedlm(io, iteration_result, ' ')
                #compute critical equilibria, and euclidean disntance to the (1, 1) equilibrium
                #given the critical parameters, compute what would be the critical equilibria using linear approximation
                #given the critical equilibria, compute what would be the critical parameters using linear approximation
                #store seed, n, d, order, pert_dir, alpha, \delta_r, \delta_r_lin, \delta_x, \delta_x_lin
            end
        end
    end
end