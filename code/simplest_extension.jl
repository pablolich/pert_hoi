using DelimitedFiles
using LinearAlgebra
include("general_perturbation_functions.jl")
# (rest of your simplest_extension.jl exactly as you provided)
#this script deals with the question: 
#how does the histogram of distances change when
#d = 2, we perturb growth rates (d_pert = 0), vary from 2 to 7, 
#the maximum perturbation magnitude is 10
#vary alpha, perturbation direction, and look at different parameter sets, indexed by their seed.

"""
Given parameters (alpha, A, B and the unperturbed equilibria x0), compute Jacobian of system.
"""
function get_jacobian(
    pars::Tuple{Vector{Float64}, Matrix{Float64}, Array{Float64,3}},
    alpha::Float64,
    x0::Vector{Float64})
    r, A, B = pars
    n = length(r)
    M = (1 - alpha) * A + alpha * sum(B[:, :, k] + B[:, k, :] for k in 1:n) .* x0'
    return M
end

"""
given parameters (vector r) and a matrix M, compute equilibria using linear approximation
"""
function get_x_star_linear(
    syst_pars::Tuple{Vector{Float64}, Matrix{Float64}, Array{Float64,3}}, 
    crit_pars::Vector{Float64}, 
    alpha::Float64, 
    x0::Vector{Float64})
    M = get_jacobian(syst_pars, alpha, x0)
    return -inv(M)*crit_pars
end

"""
given parameters (vector r) and a matrix M, compute equilibria using linear approximation
"""
function get_r_star_linear(
    syst_pars::Tuple{Vector{Float64}, Matrix{Float64}, Array{Float64,3}}, 
    crit_x::Vector{Float64}, 
    alpha::Float64, 
    x0::Vector{Float64})
    M = get_jacobian(syst_pars, alpha, x0)
    return -inv(M)*crit_x
end

"""
wrapper for readibility; runs in series all the functions necessary to prepare the system for 
feasibility boundary computation
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
system syst given variables vars for subsequent storing. 
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

nsim = 1  # Number of systems to look at
d = 2; @var α[1:d]  # Degree of polynomials to solve
pert_size = 10.0  # Maximum perturbation
d_pert = 0  # Order of parameters to perturb
n_perts = 3  # Number of perturbations
alpha_vec = [0.1, 0.9]  # Values of relative interaction strength

# Loop over n from 3 to 7
for n in 2:2
    @var x[1:n]  # Number of equations; corresponding variables
    init_sol = repeat([1.0 + 0.0im])  # Set initial solutions as two arrays
    # Build skeleton polynomial system for curren n (and d)
    ref_eqs, coeffs_mat = get_ref_polynomials(x, d, n)
    pert_dirs = points_hypersphere(n, 1.0, n_perts, false)  # Perturbation directions (try to sample regularly)
    open("../data/results_simplest_extension_n_$(n).csv", "a") do io
        # Open another file for saving all parameters for each n, seed_i
        open("../data/parameters_n_$(n).csv", "a") do param_io
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
                coeff_deg_matrix = coefficients_degree_mat(System(eqs, x), vcat(variables(syst), parameters(syst)))
                seed_column = fill(seed_i, size(coeff_deg_matrix, 1))

                println("Writing coefficients for seed $seed_i")
                # Save coefficient-degree matrix to a se    parate file named by `n` and `seed_i`
                writedlm(param_io, hcat(seed_column, coeff_deg_matrix), ' ')

                # Store the results of the simulation in the main CSV file
                for pert_i in 1:n_perts 
                    end_parameters = initial_pars .+ pert_size .* pert_dirs[pert_i,:]

                    for alpha_i in alpha_vec  # Loop through each relative strength value
                        syst_alpha = evaluate_pars(syst, α, [1-alpha_i, alpha_i])
                        println("SIMULATION: ", seed_i, " parameter set: ", pert_i, " n: ", n, " relative strength: ", alpha_i)
                        pars_crit, xstar_crit, flag = findparscrit(syst_alpha, init_sol, initial_pars, end_parameters)
                        if pars_crit == -1
                            return pars, pert_size, pert_dirs[pert_i,:], alpha_i
                        end
                        #compute euclidean distance between unperturbed and critical parameters
                        δₚ = norm(pars_crit .- initial_pars)
                        #compute euclidean distance between unperturbed and critical equilibrium
                        δₓ = norm(xstar_crit .- ones(n))  # Distance to (1, ..., 1)

                        # Linear approximations (equilibrium given critical parameters)
                        xcrit_lin = get_x_star_linear(pars, pars_crit, alpha_i, real(init_sol))
                        #(parameters given critical equilibrium)
                        rcrit_lin = get_r_star_linear(pars, real(xstar_crit), alpha_i, real(init_sol))
                        
                        # compare critical parameters and equilibrium with linear versions
                        δₓₗᵢₙ = norm(xstar_crit .- xcrit_lin)
                        δₚₗᵢₙ = norm(initial_pars .- rcrit_lin)

                        # Save results
                        iteration_result = [seed_i, n, d, pert_size, pert_i, alpha_i, δₚ, δₚₗᵢₙ, δₓ, δₓₗᵢₙ]
                        writedlm(io, iteration_result', ' ')  # Transpose to make it a row
                    end
                end
            end
        end
    end
end