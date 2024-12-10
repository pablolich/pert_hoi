#functions to perturb any parameter in a glv model of any order

using HomotopyContinuation
using IterTools
using Random

"""
generate random points on the surface of a n-dimensional hypersphere of radius rho.
when dimension is 2, the points are evenly distributed. When dimension is greater than 2 but smaller than 5
the points are loaded from the results of the python optimization that solves the Thompson problem in n dimensions
"""
function points_hypersphere(dim::Int, rho::Float64, num_points::Int, load = true)
    #no load from file, sample points at random and put them on the sphere
    if dim == 2
        points = zeros(num_points, 2)  # Initialize matrix to store points
        for i in 1:num_points
            theta = 2π * (i - 1) / num_points  # Calculate angle for current point
            points[i, 1] = rho * cos(theta)  # Calculate x coordinate
            points[i, 2] = rho * sin(theta)  # Calculate y coordinate
        end
        return points
    else
        if load == true
            return readdlm("positions_n_"*string(dim)*".csv", ',')
        else
            points = randn(num_points, dim)  # Generate random points in dim-dimensional space
            norms = [norm(points[i,:]) for i in 1:num_points]  # Calculate norms of each point
            scaled_points = rho * (points ./ norms)  # Scale points to lie on the surface of the sphere of radius rho
            return scaled_points
        end
    end
end

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

"""
given a polynomial expression and a vector its coefficients, substitute the 
"""
function symb2num(equations::Vector{Expression}, 
    coeffs::Vector{Variable},
    var_inds::Vector{Int64}, 
    var_values::Union{Vector{Real}, Vector{Complex{Float64}}})
    for i in 1:length(equations)
        equations[i] = subs(equations[i], coeffs[var_inds] => var_values)
    end
    return equations
end


"""
get a set of n dense polynomials of degree d with symbolic coefficients
"""
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
    is_feasible(r::PathResult)

return true if all components of r are positive
"""
function is_feasible(r::PathResult)
    return all(real(r.solution) .> 0) 
end

"""
    stopatnonfeasible(r::PathResult)

return true if r is not feasible
"""
function stopatnonfeasible(r::PathResult)
    #check if solution  is real
    println("solution is: ", r.solution)
    if !is_real(r) || !is_feasible(r)
        return true
    else 
        return false
    end
end

"""
given a parametrized system of polynomials, an initial set of parameters, and a final set of parameters, 
solve such system using a parameter homotopy. the solution stops the moment one of the components looses positivity
either because it becomes negative, or because it becomes complex.
"""
function solve_parametrized(syst::System, initial_pars::Vector{Float64}, final_pars::Vector{Float64}, 
                            start_solutions::Vector{Vector{ComplexF64}})
    #create the initial equilibrium
    res = solve(syst, start_solutions; 
                start_parameters = initial_pars,
                target_parameters = final_pars,
                catch_interrupt = false,
                stop_early_cb = stopatnonfeasible)
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
            slices = repeat([Colon()], j-1)
            #get corresponding tensor given order j
            T_j = (pars[j])[slices..., i]
            #create all index combinations
            index_combinations = get_inds_tensor(n, j-1)
            if j == 1
                eq_i += T_j #treat zeroth order term separately
            else
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
	    constrained_tensor[slices...,i] .= -slices_sum[i]/sum(tensor[slices...,i]) .* tensor[slices...,i]
    end
    return constrained_tensor
end

"""
sample parameters with equilibrium preserving constraints
"""
function sample_parameters(n::Int64, order::Int64, rng::MersenneTwister)
    #sample growth rates
    pars = []
    r = randn(rng, n)
    push!(pars, r)
    for i in 2:order
        push!(pars, apply_constrain(randn(rng, repeat([n], i)...), r))
    end
    pars = tuple(pars...)
    return pars
end

"""
    get_ts_and_xs(tracker_homotopy::Tracker, start_solution::AbstractVector, t₀::Float64=1.0, t₁::Float64=0.0)

given a ParameterHomotopy, an initial solution (start solution), and initial and final values of t (t0 and t1, respectively),
construct a tracker and return the values of t and x (solution of the system) as we go from t=1 to t=0
"""
function get_ts_and_xs(tracker_homotopy::Tracker, start_solution::AbstractVector, t₀::Float64=1.0, t₁::Float64=0.0)
    Xs = Vector{ComplexF64}[]
    Ts = Vector{Foat64}[]
    for (x, t) in iterator(tracker_homotopy, start_solution, t₀, t₁)
        push!(Xs, x)
        push!(Ts, t)
    end
    return Ts, Xs
end

"""
given start parameters, target parameters, and set of t and x values, construct the vector of parameters corresponding 
to the systems we traverse when doing the HomotopyContinuation
"""
function get_parameters_along_path(tvector::Vector{Float64}, initial_parameters::Vector{Float64}, target_parameters::Vector{Float64})
    #initialize parameter matrix
    ncol = length(tvec)
    nrow = length(start_parameters)
    parameter_matrix = zeros(nrow, ncol)
    for j in 1:ncol
        parameter_matrix[:,j] = tvector[i] * initial_parameters + (1 - tvector[i]) * target_parameters     
    end
    return parameter_matrix
end

"""
given the values of parameters and equilibria, determine what are the parameter values leading to negative or
complex equilibria
"""
function get_parameters_at_boundary(solution_matrix::Vector{ComplexF64}, parameter_matrix::Matrix{Float64}, tvector::Vector{Float64})
    ntsteps = length(tvector)
    critical_parameters = Vector{Float64}[]
    critical_row = 0
    for i in 1:ntsteps
        solution_i = solution_matrix[i]
        parameters_i = parameter_matrix[:, i]
        is_solution_feasible = all(x -> real(x) > 0 && imag(x) == 0, solution_i)
        if !is_solution_feasible
            return parameters_i[:, i-1]
        end
    end
end


#test all functions 

n = 2
d = 3
rng = MersenneTwister(1)

pars = sample_parameters(2, 3, rng)

#get perturbations on the unit sphere

perts = points_hypersphere(2, 1.0, 10)
rhomax = 0.5

@var  x[1:2]
@var α[1:2]

eqs = build_glvhoi(pars, x)
ref_eqs, coeffs_mat = get_ref_polynomials(x, 2, 2, coeff_name = :c)
inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
eqs_inter_str = parametrize_stengths(eqs_inter, x, α)
syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

#evaluate system at α = 0.5

#set up a particular solution to this system
start_solutions = [[1.0 + 0.0im, 1.0 + 0.0im]]
initial_parameters = vcat(r, [1., 0.])
end_parameters = vcat(r .+ rhomax*perts[8,:], [0.5, 0.5])

ct = Tracker(CoefficientHomotopy(syst; start_coefficients = initial_parameters,
                                       target_coefficients = end_parameters))

Xs = Vector{ComplexF64}[]
Ts = []
Ps = []

for (x, t, p) in iterator(ct, [-1.0], 1.0, 0.0)

push!(Xs, x)
push!(Ts, t)

end


# #create the initial equilibrium
# res = solve(syst, start_solutions;
#             start_parameters = initial_parameters,  
#             target_parameters = end_parameters,
#             catch_interrupt = false,
#             stop_early_cb = stopatnonfeasible)
# println("Solutions: ", solutions(res))

