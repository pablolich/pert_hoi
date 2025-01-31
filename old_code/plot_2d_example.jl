#code to obtain the feasibility domain as acurately as possible, in 2 dimensions, 
#and compare it to the feasibility radius obtain by looking at the circumference

using HomotopyContinuation
include("functions.jl")

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
given parameters, variables, solve the system using parameter homotopy
"""
function solve_system(x::Vector{Variable}, r::Vector{Variable}, pars::Tuple, pert_pars::Vector{Float64})
    n = length(pars[2]) #use growth rates to get number of species
    syst = System(buildglvhoi(pars, x, true), parameters = r) #build system keeping r as implicit
    res = solve(syst, [ones(n)]; 
                start_parameters = pars[2],
                target_parameters = pert_pars)
    #check number of solutions found
    n_sols = nsolutions(res)
    if n_sols == 0
        solmat = Matrix{Float64}(undef, 0, n)
    else
        solmat = mapreduce(permutedims, vcat, real_solutions(res)) 
    end
    return vec(solmat)
end

"""
function which root we want to find (evaluated by the bisection method)
"""
function get_min_xstar(x::Vector{Float64})
    if length(x) == 0
        return -1.0
    else
        return minimum(x)
    end
end

"""
Implement bisection method
"""
function bisection_method(x::Vector{Variable}, syst_pars::Vector{Variable}, pars::Tuple, direction::Vector{Float64}, x0::Float64, x1::Float64,
                          tol::Float64=1e-9)
    global xhalf = bisect(x0, x1)
    dist = abs(x0 - x1)
    while dist > tol
        #form perturbation
        r0 = pars[2]
        rnew = pars[2] .+ xhalf*direction
        #calculate equilibrium at xhalf
        xstar = solve_system(x, syst_pars, pars, rnew)
        #evaluate if the minimum is positive or negative
        value_eq = get_min_xstar(xstar)
        x0, x1 = getlims(x0, x1, xhalf, value_eq)
        dist = abs(x0 - x1)
        return bisection_method(x, syst_pars, pars, direction, x0, x1)
    end
    return xhalf
end

"""
get the actual critical perturbation boundary for which feasibility is lost
"""
function get_critical_domain(x::Vector{Variable}, r::Vector{Variable}, parameters::Tuple, n_perturbations::Int64, rperts::Matrix{Float64}, rhomax::Float64=2.0)
    #initialize boundary vectors
    r_boundary_mat = zeros(Float64, (n_perturbations, 2))
    #get growth rate perturbations
    for i in 1:n_perturbations
        rpert_i = rperts[i,:]
        rhocrit = bisection_method(x, r, parameters, rpert_i, 0.0, rhomax)
        #once rhocrit is found, build the vector r in that direction
        rcrit = rhocrit*rpert_i
        #store
        r_boundary_mat[i,:] = rcrit
    end
    return r_boundary_mat 
end

"""
for each system, get the boundaries for different alphas
"""
function many_boundaries(n_perturbations::Int64, alphavec::Vector{Float64}, nsims::Int64, x0::Float64=0.0, x1::Float64=2.0)
    #traverse t
    rng = MersenneTwister(1)
    n = 2
    tol = 1e-9
    rperts = get_angles(n_perturbations)
    @var x[1:n]
    @var r[1:n]
    for i in 1:nsims
        r0, A, B = sampleparameters(n, rng, 2)
        for alpha in alphavec
            println("Simulation: ", i, " for alpha: ", alpha)
            alpha = [alpha]
            pars = (alpha, r0, A, B)
            target = get_first_target(pars, 1e-4, rperts, n_perturbations, x, n, tol)
            r_domain = get_critical_domain(x, r, pars, n_perturbations, rperts)
            rcrit_dense = findmaxperturbation(x0, x1, get_angles(100), 100, pars, n, x, target, "follow")
            rcrit_sparse = findmaxperturbation(x0, x1, get_angles(3), 3, pars, n, x, target, "follow")
            #prepare other vectors to store
            alphaval_vec = repeat(alpha, n_perturbations)
            sim_vec = repeat([i], n_perturbations)
            tostore_domain = hcat(sim_vec, alphaval_vec, r_domain)
            tostore_rhocrit = [i alpha rcrit_dense rcrit_sparse]
            open("../data/results_2spp_boundaries.csv", :"a") do io
                writedlm(io, tostore_domain, ' ')
            end
            open("../data/results_2spp_rho_crits.csv", :"a") do io
                writedlm(io, tostore_rhocrit, ' ')
            end
        end
    end
end

n_perturbations = 1000
alphavec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
nsims = 100
many_boundaries(n_perturbations, alphavec, nsims)


