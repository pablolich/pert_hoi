
#functions to perturb any parameter in a glv model of any order

using HomotopyContinuation
using IterTools
using Random
using LinearAlgebra

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
given a polynomial expression, a vector of its symbolic coefficients, indices to evaluate, and 
values to which evaluate, perform evaluation.
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
    coeff_name::Symbol = :c
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
function parametrize_strengths(equations::Vector{Expression}, x::Vector{Variable}, alpha_vec::Vector{Variable})
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
Evaluate a set of parameters given a System of equations.
"""
function evaluate_pars(
    syst::System, 
    variables::Vector{Variable},
    values::Vector{Float64}
)
    #perform evaluation of desired parameters
    evaluated_equations = [subs(eq, variables => values) for eq in syst.expressions]
    #construct system retaining parameters that have not been evaluated
    pars_eval_syst = setdiff(parameters(syst), variables)
    return System(evaluated_equations, parameters = pars_eval_syst)
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
sample parameters with equilibrium preserving constraints, for n species, polynomial order, and seed given by rng
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
given start and end parameters, and a value of t, return the corresponding convex sum parameters
"""
function get_parameters_at_t(t::Float64, initial_parameters::Vector{Float64}, target_parameters::Vector{Float64})
    t * initial_parameters + (1 - t) * target_parameters
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
        parameter_matrix[:,j] = get_parameters_at_t(tvector[i], initial_parameters, target_parameters)
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

"""
    trackpositive!(tracker::Tracker, x::AbstractVector, [t₁=1.0, t₀=0.0]; kwargs...) -> (t_before, x_before)

Track a solution path using homotopy continuation, stopping when any component of the solution becomes negative.

This function records the last positive solution** before any variable turns negative. 
It can be useful when following solution paths where only positive (or biologically meaningful) states are valid.

# Arguments
- `tracker::Tracker`: A Tracker object that manages the path tracking.
- `x::AbstractVector`: Initial solution.
- `t₁::Float64`: Target homotopy parameter value (default `1.0`).
- `t₀::Float64`: Starting homotopy parameter value (default `0.0`).

# Keyword Arguments
- `ω::Float64`: Angle used for endgames (optional).
- `μ::Float64`: Scaling parameter for path tracking (optional).
- `extended_precision::Bool`: Whether to use extended precision (default `false`).
- `τ::Float64`: Maximum tracking time (default `Inf`).
- `keep_steps::Bool`: Whether to store intermediate steps (default `false`).
- `max_initial_step_size::Float64`: Maximum initial step size (default `Inf`).
- `debug::Bool`: Whether to output debug information during tracking (default `false`).

# Returns
- `t_before::Float64`: Homotopy parameter value just before encountering a negative component.
- `x_before::Vector{Float64}`: Solution vector just before encountering a negative component.

# Behavior
- Initializes the tracker.
- Continuously steps along the solution path.
- Monitors the solution after each step.
- If any coordinate becomes negative (real part), returns the last solution before this happens.
- If the path ends without any negative components, returns the final solution.

"""
function trackpositive!(
    tracker::Tracker,
    x::AbstractVector,
    t₁ = 1.0,
    t₀ = 0.0;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    extended_precision::Bool = false,
    τ::Float64 = Inf,
    keep_steps::Bool = false,
    max_initial_step_size::Float64 = Inf,
    debug::Bool = false,
)
    #initialize tracker
    init!(
        tracker,
        x,
        t₁,
        t₀;
        ω = ω,
        μ = μ,
        extended_precision = extended_precision,
        τ = τ,
        keep_steps = keep_steps,
        max_initial_step_size = max_initial_step_size,
    )
    #record initial t and equilibrium
    tbefore = t₁
    xbefore = x
    while is_tracking(tracker.state.code)
        #record quantities before and after step
        tbefore = tracker.state.t
        xbefore = copy(tracker.state.x)
        step!(tracker, debug)
        #after a tracking step, check if any of the components are negative.
        if any(real(tracker.state.x) .< 0)
            #return the component before negative solution was found
            return real(tbefore), real(xbefore)
        end
    end
    return real(tbefore), real(xbefore)
end

"""
    equivalent of track, but only for positive solutions
"""
@inline function trackpositive(tracker::Tracker, x, t₁ = 1.0, t₀ = 0.0; kwargs...)
    tbefore, xbefore = trackpositive!(tracker, x, t₁, t₀; kwargs...)
    TrackerResult(tracker.homotopy, tracker.state)
end

"""
Given a model syst with an equilibrium at initialsol (feasible), for parameters initial_parameters,
compute parameters for which feasibility is lost when traversing the line 
t * initial_parameters + (1-t) * target_parameters from t = 1 to t = 0. Feasibility is lost either when
the solution becomes complex (at least one entry has non-zero imaginary part) or when solution becomes
negative (at least one entry reaches 0 up to tolerance)
"""
function findparscrit(
    syst::System,
    initialsol::AbstractVector, 
    initial_parameters::Vector{Float64}, 
    target_parameters::Vector{Float64},
    tol::Float64=1e-9, 
    rec_level::Int64=1)
    #compute the euclidean distance between parameters to determine maximum step size of homotopy
    par_dist = norm(initial_parameters - target_parameters)
    max_step_size = par_dist/10 #the maximum step size should be 10 times smaller than the parameter euclidean distance
    #create a tracker to traverse parameter space from initial to target parameters
    ct = Tracker(CoefficientHomotopy(syst; start_coefficients = initial_parameters,
                                         target_coefficients = target_parameters))
    #track along path until 
        #(1) target_parameters are reached,
        #(2) negative solution is found
        #(3) complex solution is found                               
    tbefore, xbefore = trackpositive!(ct, initialsol, 1.0, 0.0, max_initial_step_size=max_step_size) #log the t and x a step before the end of the routine
    res = TrackerResult(ct.homotopy, ct.state) #form a TrackerResult
    #store some tracking output
    retcode = res.return_code #gives whether tracking succeded, failed, or sotpped
    tfinal = real(res.t) #the value of t for after tracking routine
    neg_component = minimum(real(res.solution)) #get the most negative x of solution after tracking routine
    #depending on the status after tracking, decide how to continue:
    if retcode == :success && neg_component > 0 #tracking ended succesfully at a postive solution
        #return parameters for t = 0
        #println("Tracking succesful, t = ", tfinal)
        return target_parameters
    else #tracking stopped early or failed
        if retcode == :tracking #it stopped early because negative solutions were found
            #re-run this function to get closer to the positive-negative boundary
            while abs(minimum(xbefore)) > tol #keep running until minimum is sufficiently small
                #initial solution is last positive solution of the path
                initsol = xbefore
                #initial parameters are the ones corresponding to the last postive solution along the path
                initpars = get_parameters_at_t(tbefore, initial_parameters, target_parameters)
                #end parameters are the ones correspond to the first negative solution found
                endpars = get_parameters_at_t(tfinal, initial_parameters, target_parameters)
                rec_level += 1
                if rec_level > 100
                    #the search hasn't converged
                    println("The search did not converge")
                    return target_parameters #print the parameters at the search boundary
                end
                return findparscrit(syst, initsol, initpars, endpars, tol, rec_level)
            end
            #when tolerance is reached, we are at boundary between positive-negative solutions 
            return get_parameters_at_t(tbefore, initial_parameters, target_parameters)
        else #complex solutions were found
            while maximum(abs.(imag(res.solution))) > tol
                initsol = xbefore
                initpars = get_parameters_at_t(tbefore, initial_parameters, target_parameters)
                endpars = get_parameters_at_t(tfinal, initial_parameters, target_parameters)
                rec_level += 1
                if rec_level > 100
                    println("The search did not converge (complex branch).")
                    return target_parameters #print parameters at the search boundary
                end
                return findparscrit(syst, initsol, initpars, endpars, tol, rec_level)
            end
            #end tracking since we are at boundary between positve-complex solutions
            return get_parameters_at_t(tbefore, initial_parameters, target_parameters)
        end
    end
end

"""
process output of findparscrit and generate output to be saved
1. Flag wether boundary is with negative or with complex or wether search has not converged
3. Return both parameters and equilibrium at boundary
"""
function process_output_boundary(
    )
    return pars_crit, xstar_crit, flag
end


# #test all functions 

# n = 2
# d = 3
# rng = MersenneTwister(1)

# pars = sample_parameters(2, 3, rng)

# #get perturbations on the unit sphere

# perts = points_hypersphere(2, 1.0, 10)
# rhomax = 0.5
# r = pars[1]

# @var  x[1:2]
# @var α[1:2]

# eqs = build_glvhoi(pars, x)
# ref_eqs, coeffs_mat = get_ref_polynomials(x, d, 2, coeff_name = :c)
# inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
# eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
# eqs_inter_str = parametrize_strengths(eqs_inter, x, α)
# syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

# #evaluate system at α = 0.5



# #set up a particular solution to this system
# start_solutions = [[1.0 + 0.0im, 1.0 + 0.0im]]
# initial_parameters = vcat(r, [0.5, 0.5])
# end_parameters = vcat(r .+ rhomax*perts[3,:], [0.5, 0.5])

# ct = Tracker(CoefficientHomotopy(syst; start_coefficients = initial_parameters,
#                                        target_coefficients = end_parameters))


# s = [1, 1]
# #track solution s from 1 to 0
# res = track(ct, s, 1, 0) #equivalent to solve
# #respos = trackpositive(ct, s, 1, 0) #equivalent to solve
# #need to know how to access the solution of a tracker, to check if any component is positive
# #then set the tracker.state to failure if that is the case.

# #now i have to define my own track function, which in turn calls the new track! function that 
# #stops when a component of the solution becomes negative.

# # Xs = Vector{ComplexF64}[]
# # Ts = []
# # Ps = []

# # for (x, t, p) in iterator(ct, [-1.0], 1.0, 0.0)

# # push!(Xs, x)
# # push!(Ts, t)

# # end


# #create the initial equilibrium
# println("Solutions: ", solutions(res))


#= 
This part defines modular utilities for tracking positive solutions
in homotopy continuation problems, detecting feasibility boundaries, 
and processing outputs cleanly.

Contents:
- `is_positive`: checks positivity.
- `is_real`: checks realness.
- `decide_boundary_type`: determines loss of feasibility.
- `trackpositive!` and `trackpositive`: tracking solutions while maintaining positivity.
- `track_and_monitor!`: wrapper to track and check for loss.
- `findparscrit`: main recursive boundary search.
- `process_output_boundary`: formats final outputs.
=#

# --- Basic feasibility checks ---

@inline function is_positive(x::AbstractVector, tol::Float64)
    return minimum(abs.(real(x))) < tol
end

"""
    solve_at_boundary_parameters(syst::System, t_star::Float64,
                                 initial_parameters::Vector{Float64},
                                 target_parameters::Vector{Float64})

Given a homotopy defined by `syst` and the parameters it interpolates between,
solve the system at parameters corresponding to `t_star`.

Returns: `Result` object from `solve()` with all isolated solutions.
"""
function solve_at_boundary_parameters(
    syst::System,
    t_star::Float64,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64};
    kwargs...
)
    # Interpolate parameters at t_star
    #t_star = 1
    p_star = t_star .* initial_parameters .+ (1 - t_star) .* target_parameters

    # Fix the parameters in the system
    S_fixed = fix_parameters(syst, p_star)

    # Solve the system with fixed parameters
    result = solve(S_fixed; kwargs...)

    return result
end

function decide_boundary_type(
    syst::System,
    tracker::Tracker,
    result::TrackerResult,
    tol::Float64,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64}
)
    return_code = result.return_code
    x = result.solution

    if return_code == :success
        # Tracking was successful and solution is positive
        return :boundary

    elseif return_code == :tracking
        # Tracking was interrupted manually by trackpositive!
        if any(abs.(real.(x)) .< tol)
            println("Solution at negative boundary: ", x)
            return :negative
        else
            println("Solution became negative but above tolerance: ", x)
            return :nonconverged
        end

    elseif return_code == :terminated_step_size_too_small || return_code == :terminated_max_steps
        # Tracking didn't succeed nor found negative solutions only option left is complex bifurcation.
           return :complex
    else 
        # I don't know what this could be: print
        println("Unkown boundary type: ", return_code)
        return :nonconverged
    end
end

# --- Tracking with positivity constraint ---

"""
    trackpositive!(tracker::Tracker, x::AbstractVector, [t1=1.0, t0=0.0]; kwargs...) -> (t_before, x_before)

Track a solution path, stop when any component becomes negative.
"""
function trackpositive!(
    tracker::Tracker,
    x::AbstractVector,
    t1 = 1.0,
    t0 = 0.0;
    ω::Float64 = NaN,
    μ::Float64 = NaN,
    extended_precision::Bool = false,
    τ::Float64 = Inf,
    keep_steps::Bool = false,
    max_initial_step_size::Float64 = Inf,
    debug::Bool = false,
)
    init!(
        tracker,
        x,
        t1,
        t0;
        ω = ω,
        μ = μ,
        extended_precision = extended_precision,
        τ = τ,
        keep_steps = keep_steps,
        max_initial_step_size = max_initial_step_size,
    )
    #record initial t and equilibrium
    tbefore = t1
    xbefore = x
    while is_tracking(tracker.state.code)
        #record quantities before and after step
        tbefore = tracker.state.t
        xbefore = copy(tracker.state.x)
        step!(tracker, debug)
        #after a tracking step, check if any of the components are negative.
        if any(real(tracker.state.x) .< 0)
            #return the component before negative solution was found
            return real(tbefore), real(xbefore)
        end
    end
    return real(tbefore), real(xbefore)
end

# --- Higher level tracking monitor ---

"""
    track_and_monitor!(syst, initialsol, initial_parameters, target_parameters, tol) -> (t_before, x_before, t_after, x_after, flag)

Track once from initial to target parameters, detect feasibility loss.
Returns last positive solution, solution after tracking, and boundary flag.
"""
function track_and_monitor!(
    syst::System,
    initialsol::AbstractVector,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64},
    tol::Float64,
    max_step_ratio::Float64 = 1.
)   
    par_dist = norm(initial_parameters - target_parameters)
    #the maximum step size should be 10 times smaller than the parameter euclidean distance
    max_step_size = par_dist * max_step_ratio 

    tracker = Tracker(CoefficientHomotopy(syst; start_coefficients = initial_parameters,
                                           target_coefficients = target_parameters);
                                           options = TrackerOptions(max_steps = Int(1e6))) #increased number of steps
    t_before, x_before = trackpositive!(tracker, initialsol, 1.0, 0.0, max_initial_step_size = max_step_size) #adjusted max_initial_step_size
    result = TrackerResult(tracker.homotopy, tracker.state)
    t_after = real(result.t)
    x_after = result.solution
    #store some tracking output
    flag = decide_boundary_type(syst, tracker, result, tol, initial_parameters, target_parameters)
    println("")
    println("")
    if flag == :complex
        println("AT COMPLEX BOUNDARY")
        println("t_before: ", t_after)
    end
    return t_before, x_before, t_after, x_after, flag
end

# --- Recursive critical parameter finder ---

function findparscrit(
    syst::System,
    initialsol::AbstractVector,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64},
    tol::Float64 = 1e-9,
    rec_level::Int = 1
)
    println("Recursion level: ", rec_level)
    println("Initial parameters: ", initial_parameters)
    println("Target parameters: ", target_parameters)
    t_before, x_before, t_after, x_after, flag = track_and_monitor!(
        syst, initialsol, initial_parameters, target_parameters, tol)
    #get initial and end parameters after tracking
    initpars = get_parameters_at_t(t_before, initial_parameters, target_parameters)
    endpars = get_parameters_at_t(t_after, initial_parameters, target_parameters)
    #println("Flag: ", flag)
    if flag == :negative || flag == :complex || flag == :boundary
        # No boundary was crossed, max perturbation reached
        output = process_output_boundary(endpars, x_after, flag)
        #println("")
        #println("Finished at: ", flag, " equilibrium: ", output.xstar_crit)
        return output
    elseif rec_level > 5
        #maximum recursion level reached, return result
        output = process_output_boundary(target_parameters, x_after, flag)
        return output
    else
        return findparscrit(syst, x_before, initpars, endpars, tol, rec_level + 1)
    end
end


# --- Output formatter ---
"""
    process_output_boundary(pars_crit, xstar_crit, flag)

Formats output after finding a critical point.

Returns a NamedTuple with:
- flag: Symbol (:success, :negative, :complex, or :nonconverged)
- pars_crit: Vector of critical parameters
- xstar_crit: Vector of equilibrium at boundary
"""
function process_output_boundary(
    pars_crit::Vector{Float64},
    xstar_crit::Vector{ComplexF64},
    flag::Symbol
)
    return (
        pars_crit = pars_crit,
        xstar_crit = xstar_crit,
        flag = flag
    )
end

# #test all functions 

# n = 2
# d = 2
# rng = MersenneTwister(1)

# pars = sample_parameters(2, 3, rng)

# #get perturbations on the unit sphere

# perts = points_hypersphere(2, 1.0, 10)
# rhomax = 0.5
# r = pars[1]

# @var  x[1:2]
# @var α[1:2]

# eqs = build_glvhoi(pars, x)
# ref_eqs, coeffs_mat = get_ref_polynomials(x, d, 2, coeff_name = :c)
# inds_growth_rates = get_ind_coeffs_subs(ref_eqs[1], x, "order", [0])
# eqs_inter = parametrize_interactions(eqs, ref_eqs, x, inds_growth_rates)
# eqs_inter_str = parametrize_strengths(eqs_inter, x, α)
# syst = build_parametrized_glvhoi(eqs_inter_str, x, coeffs_mat, inds_growth_rates, α)

# #evaluate system at α = 0.5

# #set up a particular solution to this system
# start_solutions = repeat([1], n)
# initial_parameters = vcat(r, [0.5, 0.5])
# end_parameters = vcat(r .+ rhomax*perts[1,:], [0.5, 0.5])

# pars_crit, xstar_crit, flag = findparscrit(syst, start_solutions, initial_parameters, end_parameters)