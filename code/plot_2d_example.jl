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
                target_parameters = get_new_pars)
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
