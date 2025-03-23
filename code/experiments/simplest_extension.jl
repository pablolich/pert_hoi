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

"""
returns a list of random coefficients (first column) and their corresponding degree (second column) of a 
system syst given variables vars. 
"""
function coefficients_degree_mat(
    syst::System,
    vars::Vector{Variable})

    n = length(syst)  # Get the number of expressions
    coeffs_list = []  # List to store coefficients for each expression
    deg_list = []     # List to store degrees for each expression

    # Loop over each expression in the system
    for i in 1:n
        # Get the exponents and coefficients for the current expression syst.expressions[i]
        e, c = exponents_coefficients(syst.expressions[i], vars)
        
        # Extract the exponents and sum them to get the degree
        exp_monomials = e[1:n, :]
        deg_monomials = vec(sum(exp_monomials, dims = 1))
        
        # Append the results to the lists
        push!(coeffs_list, c)
        push!(deg_list, deg_monomials)
    end

    # Combine the coefficients and degrees from all expressions
    coeffs_monomials_combined = vcat(coeffs_list...)
    deg_monomials_combined = vcat(deg_list...)

    result = hcat(coeffs_monomials_combined, deg_monomials_combined)
    # Remove rows where the coefficient is 1
    result_filtered = result[result[:, 1] .!= 1, :]

    # Return the matrix of coefficients and degrees
    return result_filtered
end

#set fix parameters for simulations
nsim = 1000  # Number of simulations
d = 2; @var α[1:d]  # Degree of polynomials to solve; corresponding strengths
pert_size = 10.0  # Maximum perturbation
d_pert = 0  # Order of parameters to perturb
n_perts = 1000  # Number of perturbations
alpha_vec = [0.1, 0.9]  # Values of relative interaction strength

# Loop over n from 3 to 7
for n in 3:7
    @var x[1:n]  # Number of equations; corresponding variables
    init_sol = repeat([1], n)  # Set initial solution to 1
    # Build skeleton polynomial system for n, d
    ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n)
    pert_dirs = points_hypersphere(n, 1.0, n_perts, false)  # Perturbation directions
    open("../../data/results_simplest_extension_n_$(n).csv", "a") do io
        # Open another file for saving all parameters for each n, seed_i
        open("../../data/parameters_n_$(n).csv", "a") do param_io
            # Run all simulations
            for seed_i in 1:nsim
                rng = MersenneTwister(seed_i)
                # Sample parameters. Note that they are deterministic given n, order and rng.
                pars = sample_parameters(n, d+1, rng)
                initial_pars = pars[1]
                # Build numerical glvhoi system
                eqs = build_glvhoi(pars, x)

                # Build system parametrizing strengths and parameters to be perturbed
                syst = get_parametrized_system(eqs, ref_eqs, coeffs_mat, d_pert, x, α)

                # Calculate the coefficients and degrees for the system
                coeff_deg_matrix = coefficients_degree_mat(syst, vcat(variables(syst), parameters(syst)))
                seed_column = fill(seed_i, size(coeff_deg_matrix, 1))

                
                # Save coefficient-degree matrix to a se    parate file named by `n` and `seed_i`
                writedlm(param_io, hcat(seed_column, coeff_deg_matrix), ' ')

                # Store the results of the simulation in the main CSV file
                for pert_i in 1:n_perts 
                    end_parameters = initial_pars .+ pert_size .* pert_dirs[pert_i,:]

                    for alpha_i in alpha_vec  # Loop through each relative strength value
                        syst_alpha = evaluate_pars(syst, α, [1-alpha_i, alpha_i])
                        println("Simulation: ", seed_i, " parameter set: ", pert_i, " n: ", n, " relative strength: ", alpha_i)
                        pars_crit = findparscrit(syst_alpha, init_sol, initial_pars, end_parameters)
                        if pars_crit == -1
                            return pars, pert_size, pert_dirs[pert_i,:], alpha_i
                        end
                        δₚ = norm(pars_crit .- initial_pars)
                        
                        iteration_result = [seed_i n alpha_i pert_i δₚ]
                        writedlm(io, iteration_result, ' ')  # Store simulation results
                        # Compute critical equilibria, and Euclidean distance to the (1, 1) equilibrium
                        # Given the critical parameters, compute what would be the critical equilibria using linear approximation
                        # Given the critical equilibria, compute what would be the critical parameters using linear approximation
                        # Store seed, n, d, order, pert_dir, alpha, \delta_r, \delta_r_lin, \delta_x, \delta_x_lin
                    end
                end
            end
        end
    end
end
