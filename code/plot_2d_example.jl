#sample parameters n = 2
#pick one alpha
#get the real feasibility domain
#get critical radius with many points (1e4)
#get critical radius with a fraction of points (1e2)
#plot the r's coming from the real feasibility domain
#plot a circle with the radius equal the one obtained for 1e4 (i.e., inf)
#plot a the points used to find critical radius of 1e2 points
#if desired, plot also the x corresponding to this cases

"""
find middle point between two points
"""
function bisect(x0::Float64, x1::Float64)
    return 1/2*(x0 + x1)
end

"""
get appropriate limits based on sign of the middle point
"""
function getlims(x0::Float64, x1::Float64, xb::Float64, xmin::Float64)
    if xmin > 0
        x0 = xb
    else
        x1 = xb
    end
    return x0, x1
end

"""
get num_points equispaced in a circle of radius 1
"""
function get_angles(num_points::Int64)
    points = zeros(num_points, 2)  # Initialize matrix to store points
    for i in 1:num_points
        theta = 2π * (i - 1) / num_points  # Calculate angle for current point
        points[i, 1] = cos(theta)  # Calculate x coordinate
        points[i, 2] = sin(theta)  # Calculate y coordinate
    end
    return points
end

"""
build glv model
"""
function buildglvhoi(pars::Tuple, x::Vector{Variable}, symb_pars::Bool=false)
    n = length(pars[2])
    if symb_pars == false
        #unpack parameters
        alpha, r, A, B = pars
    else
        #unpack parameters keeping r as a variable
        alpha, _, A, B = pars
        @var r[1:n]
    end
    eqs = r + (1 .- alpha) .* A*x
    #add HOIs
    for i in 1:n
        eqs[i] += (alpha .* ( x'*B[:,:,i]*x ))[1]
    end
    return diagm(x) * eqs
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
    for i in 1:length(coeff_inds)
        num_polynomial = subs(num_polynomial, coeffs_num[i] => coeffs_ref[i])
    end
    return num_polynomial
end

"""
build new tuple of parameters given the old and the perturbation
"""
function get_new_pars(old_pars::Tuple, perturbation::Tuple)
    #get length of the tuple
    n_par_group = length(old_pars)
    for i in 1:n_par_group
        old_pars[i] += perturbation[i]
    end
    return old_pars
end

"""
given parameters, variables, solve the system using parameter homotopy
"""
function solve_system(vars, pars, pert_pars)
    n = length(pars[2]) #use growth rates to get number of species
    syst = System(buildglvhoi(pars, x, true), parameters = r) #build system keeping r as implicit
    start_solutions = [ones(n)]
    new_pars = get_new_pars(pars, pert_pars)
    res = solve(syst, start_solutions=[ones(n)]; 
                start_parameters = pars[2],
                target_parameters = new_pars)
    return
end

"""
get the actual critical perturbation boundary for which feasibility is lost
"""
function get_critical_domain(parameters::Tuple, )
    #pick a direction and increase r until feasibility is lost
    #if feasibility is not lost, then pick another (random) direction and do the same thing
    #to see determine when is feasibility lost, use bisection method.
    #in particular, find the r that makes at least one component of x 0 or complex
    #once that point is found, rotate the r slightly and perform the search again.
    #the limits of the search should be a 15% increase on both directions from the previous search
    return 
end
