#functions to perturb any parameter in a glv model of any order

using HomotopyContinuation
using IterTools
using Random

"""
given a polynomial and a mode:
if mode == "order" then return the indices of monomials order focal_monomial
if mode == "species" return the indieces of the monomials containing the variable in focal_monomial
if mode == "both" then return the intersection of the two, so only the coefficients of order focal_monomial[1], 
in which the species focal_monomial[2] is involved.
IMPORTANT: if mode is both, then the length of focal_monomial must be 2, in any other case it must be 1
"""
function get_ind_coeffs_subs(polynomial::Expression, x::Vector{Variable}, mode::String, focal_monomial::Vector{Int64})
    inds_order = []
    inds_species = []
    #get the exponents of each variable in each monomial and the associated coefficients
    exponents, coeffs = exponents_coefficients(polynomial, x)
    #summ the exponents to determine the order of each monomial
    order_monomials = vec(sum(exponents, dims = 1)) #vectorize too
    if mode == "order"  
        #find positions of monomials of order focal_monomial
        inds_order = findall(x -> x == focal_monomial[1], order_monomials) 
        return inds_order
    elseif mode == "species"
        #focal_monomial should be interpreted as an integer in [0, n], where n is the number of variables
        exps_focal_species = exponents[focal_monomial[1],:]
        #find positions where variable focal_monomial appear
        inds_species = findall(x -> x != 0, exps_focal_species)
        return inds_species
    elseif mode == "both"
        inds_order = findall(x -> x == focal_monomial[1], order_monomials)
        #focal_monomial should be interpreted as an integer in [0, n], where n is the number of variables
        exps_focal_species = exponents[focal_monomial[2],:]
        #find positions where variable focal_monomial appear
        inds_species = findall(x -> x != 0, exps_focal_species)
        return intersect(inds_order, inds_species)
    else
        throw(ErrorException("Not a valid mode to select coefficients"))
    end
end

"""
given ref_polynomial (a reference polynomial) with symbolic coefficients,  num_polynomial (a polynomial with 
numeric coefficients), and a vector of coefficient indices: transform the numerical coefficients of num_polynomial 
into the symbolic coefficients of ref_polynomial
"""
function num2symb(ref_polynomial::Expression, num_polynomial::Expression, x::Vector{Variable}, coeff_inds::Vector{Int64})
    #puts the coefficients of ref_polynomial and num_polynomial into vectors
    coeffs_ref = coefficients(ref_polynomial, x)
    coeffs_num = coefficients(num_polynomial, x)
    #substitute the desired coefficients into the numerical polynomial
    for i in coeff_inds
        num_polynomial = subs(num_polynomial, coeffs_num[i] => coeffs_ref[i])
    end
    return num_polynomial
end

function get_ref_polynomials(
    vars::AbstractVector{<:Union{Variable,Expression}},
    d::Integer,
    n::Integer;
    homogeneous::Bool = false,
    coeff_name::Symbol = gensym(:c),
)
    eqs = Vector{Expression}()
    M = monomials(vars, d; affine = !homogeneous)
    cmat = reshape([Variable(coeff_name, i, j) for i in 1:n for j in 1:length(M)], length(M), n)
    for i in 1:n
        c_i = cmat[:,i]
        append!(eqs, sum(c_i .* M))
    end
    eqs, cmat
end

# function get_interaction_variables(ref_polynomials::Vector{Expression}, x::Vector{Variable}, coeff_inds::Vector{Int64})
#     n_equations = length(ref_polynomials)
#     #initialize
#     interaction_variables = Vector{Variable}()
#     for i in 1:n_equations
#         syst_parameters = coefficients(ref_polynomials[i], x)[coeff_inds]
#         variables_eq_i = [Variable(string(syst_parameters[i])) for i in 1:length(syst_parameters)]
#         push!(interaction_variables, variables_eq_i)
#     end
#     return interaction_variables
# end

"""
given a list of polynomial equations with numerical coefficients, a reference polynomial, and a vector of coefficient indices
to perturb, compute another set of polynomial equations where the numerical coefficients are now symbolic parameters
"""
function parametrize_interactions(equations::Vector{Expression}, ref_polynomials::Vector{Expression}, x::Vector{Variable}, coeff_inds::Vector{Int64})
    n_equations = length(equations)
    #initialize parametrized equations
    parametrized_equations = Vector{Expression}()
    for i in 1:n_equations
        append!(parametrized_equations, num2symb(ref_polynomials[i], equations[i], x, coeff_inds))
    end    
    return parametrized_equations
end

"""
given a set of polynomials, multiply all the monomials of a given order by the corresponding strength of that order given by
alpha_vec
"""
function parametrize_stengths(equations::Vector{Expression}, x::Vector{Variable}, alpha_vec::Vector{Variable})
    n_equations = length(equations)
    parametrized_equations = Vector{Expression}()
    for i in 1:n_equations
        #get equation of species i
        eq_i = equations[i]
        #get exponents of each monomial
        exponents, coeffs = exponents_coefficients(eq_i, x)
        #get number of coefficients corresponding to each order
        monomial_orders = vec(sum(exponents, dims = 1))
        #remove 0 because the alphas only affect the non-linear target_parameters
        monomial_orders = filter(x -> x != 0, monomial_orders)
        #get alphas corresponding to each coefficient
        alpha_extended = vcat(alpha_vec[monomial_orders], 1) #put back the 1 to fix dimensions before multiplication
        new_coeffs = alpha_extended .* coeffs
        new_equation = sum(monomials(x, maximum(monomial_orders)).*new_coeffs)
        append!(parametrized_equations, new_equation)
    end
    return parametrized_equations
end

"""
form a System of polynomials given a set of equations, and parameters (both interactions and strengths)
"""
function build_parametrized_glvhoi(equations::Vector{Expression}, 
                                   x::Vector{Variable},
                                   interaction_coeffs::Matrix{Variable},
                                   interactions_to_perturb::Vector{Int64},
                                   alpha_vec::Vector{Variable})
    #vectorize interaction variables to be perturbed
    interaction_parameters = vec(interaction_coeffs[interactions_to_perturb,:])
    println("parameters: ", vcat(interaction_parameters, alpha_vec))
    #build parametrized system
    return System(equations, 
                  variables = x,
                  parameters = vcat(interaction_parameters, alpha_vec))
end

"""
given a parametrized system of polynomials, an initial set of parameters, and a final set of parameters, 
solve such system using a parameter homotopy
"""
function solve_parametrized(syst::System, initial_pars::Vector{Float64}, final_pars::Vector{Float64}, 
                            start_solutions::Vector{Vector{ComplexF64}})
    #create the initial equilibrium
    res = solve(syst, start_solutions; 
                start_parameters = initial_pars,
                target_parameters = final_pars)
    #return all roots
    solvecs = solutions(res)
    if isempty(solvecs) #solution became complex, but parameter homotopy could not follow it
        solmat = Matrix{ComplexF64}(undef, 0, n)
    else
        solmat = mapreduce(permutedims, vcat, solvecs)
    end
    return solmat
end

"""
given a vector nxnx...xn (d times) of dimensions of a tensor, generate all the index positions of its elements
"""
function get_inds_tensor(n::Int64, d::Int64)
    dimensions = repeat([n], d)
    ranges = [1:d for d in dimensions]
    indices = collect(IterTools.product(ranges...))
    return vec(indices)
end

"""
given a set of parameters in the matrix space of (r, A, B, ...) construct a the glv polynomial equations
"""
function build_glvhoi(pars, x)
    n = length(x)
    #get order of glvhoi model
    d = length(pars) #subtract the 0th order term
    eqs = Vector{Expression}()
    for i in 1:n
        #a way to recursively increment the order
        #a way to set the alphas
        eq_i = 0.0
        #loop through orders
        for j in 1:d
            #treat order 0 and 1 separately
            if j == 1
                eq_i += (pars[1])[i]
            elseif j == 2
                eq_i += (pars[2]*x)[i]
            #for order 2 and above, do it in general
            else
                slices = repeat([Colon()], j-1)
                #get corresponding tensor given order j
                T_j = (pars[j])[slices..., i]
                #create all index combinations
                index_combinations = get_inds_tensor(n, j-1)
                for ind_comb in 1:length(index_combinations)
                    index_comb = index_combinations[ind_comb]
                    vec_vars = [x[i] for i in index_comb]
                    prod_vars = prod(vec_vars)
                    tensor_element = T_j[index_comb...]
                    eq_i += tensor_element*prod_vars
                end
            end
        end
        push!(eqs, eq_i)
    end
    return eqs
end

"""
given a tensor, apply the constrain in slices_sum
"""
function apply_constrain(tensor::Array{Float64}, slices_sum::Vector{Float64})
    #make the slices
    #add up all the slices except the ith, and divide by that
    #test by actually summing all the slices and checking that it adds up to the desired
    #stronger test if equilibrium of glvhois contains the vector of ones.
    S = size(tensor)
    d = length(S)
    n = S[1]
    slices = repeat([Colon()], d-1)
    constrained_tensor = zeros(repeat([n], d)...)
    for i in 1:n
	    constrained_tensor[slices...,i] .= -slices_sum/sum(tensor[slices...,i]) .* tensor[slices...,i]
    end
    return constrained_tensor
end

"""
sample parameters with equilibrium preserving constraints
"""
function sample_parameters(n::Int64, order::Int64, rng::MersenneTwister)
    #sample growth rates
    r = randn(rng, n)
    pars = (r)
    for i in 1:order
        pars = tuple(pars..., apply_constrain(randn(rng, repeat([n], order)...), r))
    end
    return pars
end


#test all functions 

n = 2
r = randn(2)
A = randn(2,2)
B = randn(2,2,2)

pars = (r, A, B)

@var  x[1:2]
@var α[1:2]

eqs = build_glvhoi(pars, x)
ref_eqs, coeffs_mat = get_ref_polynomials(x, 2, 2, coeff_name = :c)
eqs_inter = parametrize_interactions(eqs, ref_eqs, x, [5,6])
eqs_inter_str = parametrize_stengths(eqs_inter, x, α)
syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, [5, 6], α)


########################################################################################################
#FUNCTIONS I SCAVANGED FROM OLD CODE THAT MIGHT BE USEFUL here
########################################################################################################
"""
    is_feasible(r::PathResult)

return true if all components of r are positive
"""
function is_feasible(r::PathResult)
    return all(real(r.solution) .> 0) 
end

"""
    stopatreal(r::PathResult)

    return true if r is real
"""
function stopatreal(r::PathResult)
    #check if solution  is real
    if is_real(r)
        return true
    else
        return false
    end
end

"""
    stopatfeasible(r::PathResult)

return true if r is feasible
"""
function stopatfeasible(r::PathResult)
    #check if solution  is real
    if is_real(r)
        #check if its feasible
        if is_feasible(r)
            return true
        else
            return false
        end
    else
        return false
    end
end
