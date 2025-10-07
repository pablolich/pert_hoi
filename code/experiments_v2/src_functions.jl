# ─────────────────────────────────────────────────────────────────────────────
# src_functions.jl
# ─────────────────────────────────────────────────────────────────────────────

module U

using Distributions, Random, LinearAlgebra, Serialization
using Base.Iterators: product


using HomotopyContinuation


export sample_parameters_dirichlet,
       alpha_k_ratio,
       baseline_robustness_GLV,
       points_hypersphere,
       build_glvhoi,
       get_parametrized_system,
       findparscrit,
       plug_parameters_into_system,
       pairwise_distances,
       append_pairwise_distances!,
       append_no_pairs!,
       build_square_mesh


# Constrain A so that sum_j A_ij = -r_i
function constrain_matrix_dirichlet(n::Int, r::Vector, rng::AbstractRNG)
    A = zeros(n, n)
    for i in 1:n
        w = rand(rng, Dirichlet(fill(1.0, n)))  # uniform Dirichlet over n entries
        A[:, i] .= -r[i] * w
    end
    return A
end

# Constrain B so that sum_{j,k} B_{ijk} * x*_j * x*_k = -r_i
function apply_constrain_dirichlet(n::Int, order::Int, xstar::Vector, r::Vector, rng::AbstractRNG)
    T = zeros(Float64, ntuple(_ -> n, order)...)
    inds = CartesianIndices(ntuple(_ -> n, order - 1))
    for i in 1:n
        w = rand(rng, Dirichlet(fill(1.0, length(inds))))
        for (j, idx) in enumerate(inds)
            prod_x = prod(xstar[k] for k in Tuple(idx))
            T[Tuple(idx)..., i] = -r[i] * w[j] / prod_x
        end
    end
    return T
end

# Full parameter sampling: (r, A, B₂,…,B_order)
function sample_parameters_dirichlet(n::Int, order::Int, rng::AbstractRNG)
    xstar = ones(n)
    r = randn(rng, n)                       # growth rates
    A = constrain_matrix_dirichlet(n, r, rng)

    pars = Any[r, A]
    for o in 3:order
        B = apply_constrain_dirichlet(n, o, xstar, r, rng)
        push!(pars, B)
    end
    return pars
end
 
# Constraint matrix builder
"""
    build_C(ℓ, m, n)

Construct the 2n×D matrix C enforcing:
 1) For each i=1:n, sum over A[i, j] plus r[i] == 0
 2) For each i=1:n, sum over B[:,:,i] plus r[i] == 0
Where D = n + m*n + ℓ*m*n.
"""
function build_C(ℓ::Int, m::Int, n::Int)
    D = n + m*n + ℓ*m*n
    C = zeros(Float64, 2n, D)

    # 1) sum_j A[i,j] + r[i] == 0  (row sums of A)
    for i in 1:n
        C[i, i] = 1.0
        for j in 1:n
            col = n + (j-1)*m + i
            C[i, col] = 1.0
        end
    end

    # 2) sum_{j,k} B[j,k,i] + r[i] == 0
    offsetB = n + m*n
    for i in 1:n
        row = n + i
        C[row, i] = 1.0
        for k in 1:m, j in 1:ℓ
            lin = j + (k-1)*ℓ + (i-1)*ℓ*m
            C[row, offsetB + lin] = 1.0
        end
    end

    return C
end

# Modified hit-and-run sampler using a truncated Normal(0,1) step
"""
    hit_and_run_sample(C, L, U, x0, ℓ, m, n)

Perform one hit-and-run move sampling λ from a truncated Normal(0,1) on the feasible chord.
Returns (r, A, B)
"""
function hit_and_run_sample(
    C::AbstractMatrix{<:Float64},
    L::AbstractVector{<:Float64},
    U::AbstractVector{<:Float64},
    x0::AbstractVector{<:Float64},
    ℓ::Int, m::Int, n::Int;
    max_tries::Int = 100
)
    # Precompute nullspace basis
    F     = svd(C; full=true)
    tol   = maximum(size(C)) * eps(eltype(C)) * maximum(F.S)
    rankC = count(>(tol), F.S)
    Q     = F.V[:, rankC+1:end]
    D     = length(x0)

    for attempt in 1:max_tries
        # random direction in ker(C)
        u = Q * randn(size(Q,2))
        nu = norm(u)
        if nu < 1e-12
            continue
        end
        u ./= nu

        # find feasible interval
        λ_min, λ_max = -Inf, Inf
        for j in 1:D
            if abs(u[j]) > 1e-12
                a = (L[j] - x0[j]) / u[j]
                b = (U[j] - x0[j]) / u[j]
                lo, hi = min(a,b), max(a,b)
                λ_min = max(λ_min, lo)
                λ_max = min(λ_max, hi)
            end
        end
        if λ_max < λ_min || !isfinite(λ_min) || !isfinite(λ_max)
            continue
        end

        # sample λ from truncated Normal(0,1)
        λ = rand() * (λ_max - λ_min) + λ_min
        x = x0 .+ λ .* u

        # split back into (r, A, B)
        r =     x[1                 : n]
        A = reshape(x[n+1             : n + m*n], m, n)
        B = reshape(x[n + m*n + 1     : end],       ℓ, m, n)

        return r, transpose(A), B
    end
    error("Failed to find a valid hit-and-run step after \$max_tries tries.")
end

"""
    generate_and_save_parameter_sets(
        nsim::Int, ℓ::Int, m::Int, n::Int, folder::String;
        lower, upper, prefix
    )

Run `nsim` hit-and-run samples and save each tuple (r, A, B) in its own file.
Files are named `prefix_seed#.bin` and contain the serialized `(r, A, B)` tuple.
"""
function generate_and_save_parameter_sets(
    nsim::Int, ℓ::Int, m::Int, n::Int, folder::String;
    lower::Float64 = -10.0,
    upper::Float64 = 10.0,
    prefix::String = "param"
)
    # build shared problem data
    D = n + m*n + ℓ*m*n
    C = build_C(ℓ, m, n)
    L = fill(lower, D)
    U = fill(upper, D)
    x0 = vcat(zeros(n), zeros(m*n), zeros(ℓ*m*n))

    # ensure folder exists
    isdir(folder) || mkpath(folder)

    for seed in 1:nsim
        seed =3
        rng = MersenneTwister(seed)
        # single sample
        r, A, B = hit_and_run_sample(C, L, U, x0, ℓ, m, n)
        # save tuple
        filepath = joinpath(folder, "$(prefix)_seed$(seed).bin")
        open(filepath, "w") do io
            serialize(io, (r, A, B))
        end
    end
end

"""
create a 2D square mesh around 0 of size boxsize containing sqrt(npoints)
"""
function build_square_mesh(npoints::Int, boxsize::Float64)
    lengthvec = Int(round(sqrt(npoints)))
    outer_vec = collect(range(-boxsize/2,boxsize/2, length = lengthvec))
    return [collect(Iterators.product(outer_vec, outer_vec))...]
end

# — α so that Var[(1−α)A] = k·Var[αB]
alpha_k_ratio(σA::Float64, σB::Float64, k::Float64=1.0) = 1/(1 + sqrt(k)*σB/σA)

# — analytic‐GLV baseline: x* = -A⁻¹r ⇒ feasibility lost when any coord hits zero
function baseline_robustness_GLV(A::AbstractMatrix, r::AbstractVector,
                                 pert_dirs::AbstractMatrix, Δmax::Float64)
    invA = inv(A)                    # O(n³) once
    x0   = -invA * r                 # baseline equilibrium (here all ones)
    tmax = Float64[]                 # output vector, length m

    for v in eachrow(pert_dirs)
    # how the equilibrium shifts per unit Δ:
    y = invA * v                   # so x*(Δ) = x0 .- Δ*y

    # only species with y[i]>0 will decline as Δ grows
    # find for each such i the root Δ_i = x0[i]/y[i]
    ts = [ x0[i] / y[i] for i in eachindex(r) if y[i] > 0 ]

    # record the first extinction Δ (or Δmax if none)
    push!(tmax, isempty(ts) ? Δmax : minimum(ts))
    end

    return tmax

end

"""
generate random points on the surface of a n-dimensional hypersphere of radius rho.
when dimension is 2, the points are evenly distributed. When dimension is greater than 2 but smaller than 5
the points are loaded from the results of the python optimization that solves the Thompson problem in n dimensions
"""
function points_hypersphere(dim::Int, rho::Float64, num_points::Int, load = true, rng::AbstractRNG = MersenneTwister(42))
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
given a vector nxnx...xn (d times) of dimensions of a tensor, generate all the index positions of its elements
"""
function get_inds_tensor(n::Int64, d::Int64)
    dimensions = repeat([n], d)
    ranges = [1:n for _ in 1:d]
    indices = collect(product(ranges...))
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
    println("order_monomials: ", order_monomials)
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
given start and end parameters, and a value of t, return the corresponding convex sum parameters
"""
function get_parameters_at_t(t::Float64, initial_parameters::Vector{Float64}, target_parameters::Vector{Float64})
    t * initial_parameters + (1 - t) * target_parameters
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
            #println("Solution at negative boundary: ", x)
            return :negative
        else
            #println("Solution became negative but above tolerance: ", x)
            return :nonconverged
        end

    elseif return_code == :terminated_step_size_too_small || return_code == :terminated_max_steps
        # Tracking didn't succeed nor found negative solutions only option left is complex bifurcation.
           return :complex
    else 
        # I don't know what this could be: print
        #println("Unkown boundary type: ", return_code)
        return :unknown
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
    if flag == :complex
        # println("AT COMPLEX BOUNDARY")
        # println("t_before: ", t_after)
    end
    return t_before, x_before, t_after, x_after, flag
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

# --- Recursive critical parameter finder ---

function findparscrit(
    syst::System,
    initialsol::AbstractVector,
    initial_parameters::Vector{Float64},
    target_parameters::Vector{Float64},
    tol::Float64 = 1e-9,
    rec_level::Int = 1
)
    #println("Recursion level: ", rec_level)
    #println("Initial parameters: ", initial_parameters)
    #println("Target parameters: ", target_parameters)
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


# --- Parameter substitution ---
"""
    plug_parameters_into_system(syst::System, pvals::AbstractVector{<:Real}) -> System

Return a parameter-free `System` by plugging `pvals` for all parameters of `syst`.
"""
function plug_parameters_into_system(syst::System, pvals::AbstractVector{<:Real})
    pars = parameters(syst)
    @assert length(pars) == length(pvals) "Expected $(length(pars)) parameter values, got $(length(pvals))."

    # 1) Evaluate the equations at the given parameter values
    eqs_num = evaluate(expressions(syst), pars => pvals)  # Vector{Expression}

    # 2) Rebuild a System with the same variables, but now no free parameters
    sys_num = System(eqs_num; variables = variables(syst))

    # (optional) sanity check: should be parameter-free now
    @assert isempty(parameters(sys_num)) "Unexpected remaining parameters after substitution."

    return sys_num
end

# --- Pairwise Euclidean distances ---
"""
    pairwise_distances(xcrit, sols; atol=1e-6)
        -> Vector{Tuple{Int,Int,Float64,String,String}}

Compute distances only from the critical solution `xcrit` (i = 0, always "real")
to each solution `sols[j]` returned by `solutions(solve(...))`.

- Distance uses only the real parts.
- j is flagged "real" if max(abs.(imag.(sols[j]))) ≤ atol, otherwise "complex".
- Returns tuples: (sol_i=0, sol_j=j, dist, sol_i_type="real", sol_j_type).
"""
function pairwise_distances(
    xcrit::AbstractVector,
    sols::AbstractVector{<:AbstractVector{<:Complex}};
    atol::Float64 = 1e-6,
)
    xr = real.(xcrit)
    out = Vector{Tuple{Int,Int,Float64,String,String}}()
    @inbounds for (j, sj) in pairs(sols)
        y  = real.(sj)
        dj = norm(xr .- y)
        jtype = maximum(abs.(imag.(sj))) <= atol ? "real" : "complex"
        push!(out, (0, j, dj, "real", jtype))
    end
    return out
end

"""
    append_pairwise_distances!(path,
                               seed, pert_i, αval,
                               pairs; nvars=nothing, nsols=nothing)

Append rows with pairwise distances to `path`. Columns:
seed,pert_i,alpha,sol_i,sol_j,dist,sol_i_type,sol_j_type,nvars,nsols
"""
function append_pairwise_distances!(path::AbstractString,
                                    seed::Int, pert_i::Int, αval::Real,
                                    pairs::Vector{Tuple{Int,Int,Float64,String,String}};
                                    nvars::Union{Nothing,Int}=nothing,
                                    nsols::Union{Nothing,Int}=nothing)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
    newfile = !isfile(path)
    open(path, newfile ? "w" : "a") do io
        if newfile
            println(io, "seed,pert_i,alpha,sol_i,sol_j,dist,sol_i_type,sol_j_type,nvars,nsols")
        end
        for (i, j, d, itype, jtype) in pairs
            println(io, "$(seed),$(pert_i),$(Float64(αval)),$(i),$(j),$(d),$(itype),$(jtype),",
                        isnothing(nvars) ? "" : string(nvars), ",",
                        isnothing(nsols) ? "" : string(nsols))
        end
    end
    return nothing
end

"""
    append_no_pairs!(path, seed, pert_i, αval; nvars=nothing, nsols=nothing)

Append a single row when there are no distances to record
(e.g., no solutions found). Columns match the pairwise file.
"""
function append_no_pairs!(path::AbstractString, seed::Int, pert_i::Int, αval::Real;
                          nvars::Union{Nothing,Int}=nothing, nsols::Union{Nothing,Int}=nothing)
    dir = dirname(path)
    isdir(dir) || mkpath(dir)
    newfile = !isfile(path)
    open(path, newfile ? "w" : "a") do io
        if newfile
            println(io, "seed,pert_i,alpha,sol_i,sol_j,dist,sol_i_type,sol_j_type,nvars,nsols")
        end
        # sol_i = 0 (critical), sol_i_type = real; no j
        println(io, "$(seed),$(pert_i),$(Float64(αval)),0,,,real,,",
                    isnothing(nvars) ? "" : string(nvars), ",",
                    isnothing(nsols) ? "" : string(nsols))
    end
    return nothing
end

end # module
