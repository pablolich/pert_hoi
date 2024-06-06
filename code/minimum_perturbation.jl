#here I write code to find what is the minimum perturbation (in the 2-norm sense) that turns 
#an initially feasible equilibria unfeasible.
#I start with an exhaustive exploration of the r's, to map the feasibility domain, then pick the shortest vector
#leading to extinction.
#set n=3
#Pick random parameters r0, A, B leading to a feasible equilibrium.

using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Random

include("experiment.jl")

"""
sample parameters such that ones(n) is a zero of the system
"""
function sampleparameters(n, rng)
    r0 = randn(rng)
    r = repeat([r0], n)
    A = randn(rng, (n,n))
    B = randn(rng, (n,n,n))
    Aconstr = constrainA(A, r0) #equilibrium preserving constraints
	#Bconstr = constrainB(B, r0, n) #only equilibrium preserving constraints.
    Bconstrstoch = getconstantsumB(B, n, r0) #equilibrium and tractability constraints.
    return r, Aconstr, Bconstrstoch
end

"""
transform spherical coordinates to cartesian coordinates
"""
function sph2cart(theta, phi, rho)
    x = rho*sin(phi)*cos(theta)
    y = rho*sin(phi)*sin(theta)
    z = rho*cos(theta)
    r = [x, y, z]
    return r
end

"""
transform polar coordinates to cartesian coordinates
"""
function pol2cart(theta, rho)
    x = rho*cos(theta)
    y = rho*sin(theta)
    r = [x, y]
    return r
end


"""
build glv model
"""
function buildglvhoi(pars, x)
    #unpack parameters
    alpha, r, A, B = pars
    n = length(r)
    eqs = r + (1 .- alpha) .* A*x
    #add HOIs
    for i in 1:n
	eqs[i] += (alpha .* ( x'*B[:,:,i]*x ))[1]
    end
    return diagm(x) * eqs
end

"""
function whose zero we want to find
"""
function xstarpert(r, pars, x)
    alpha, r, A, B = pars
    #get parameters
    rho = r + r
    #form system
    syst = buildglvhoi(pars, x)
    res = HomotopyContinuation.solve(syst ./ x, compile = false)
	solmat = real_solutions(res)
    ind_closest2zero = findmin(abs.(xstar))[2]
    return xstar[ind_closest2zero]
end

"""
perform minimization to find the minimum r that leads to extinction
"""
function find_r_mag(func, )
    find_zero(xstarpert, (0, 10))
    return result
end

"""
check if any row of given matrix has all posiutive components
"""
function anyrowfeasible(matrix::Matrix{T}) where T <: Real
    for row in eachrow(matrix)
        if all(x -> x > 0, row)
            return true
        end
    end
    return false
end

"""
get indices where there is a feasible equilibria
"""
function getfeasrowinds(matrix::Matrix{T}) where T <: Real
    rowind = 0
    rowinds = []
    for row in eachrow(matrix)
        rowind += 1
        if all(x -> x > 0, row)
            push!(rowinds, rowind)
        end
    end
    return rowinds
end

"""
get matrix with real solutions of system
"""
function makeandsolve(x, pars)
    #make and solve new system
    syst = System(buildglvhoi(pars, x))
    res = solve(syst ./ x, compile = false)
    solvecs = real_solutions(res)
    if nreal(res) == 0
        solmat = Array{Float64}(undef, 0, 0)
    else
        solmat = mapreduce(permutedims, vcat, solvecs)
    end
    return solmat
end

"""
identify the correct equilibrium after a small change to parameters
"""
function identifyequilibrium(eq_og, eq_candidates)
    distance_vec = []
    #get distances from candidates to the original equilibria
    for row in eachrow(eq_candidates)
        distance = norm(eq_og .- row)
        push!(distance_vec, distance)
    end
    #find minimum distance
    ind_min = argmin(distance_vec)
    return ind_min
end

"""
find parameters leading to feasibility
"""
function getfeasiblepars(n, alpha, rng)
    feasible = false
    r = []
    A = []
    B = []
    while !feasible
        r = randn(rng, n)
        A = randn(rng, (n,n))/n
        B = randn(rng, (n,n,n))/n^2
        pars = (alpha, r, A, B)
        solmat = makeandsolve(x, pars)
        feasible = anyrowfeasible(solmat)
    end
    return r, A, B
end

"""
given two vectors, calculate area of subtended triangle
"""
function getareatriangle(r0, r1)
    #transform to length 3
    r0 = vcat(r0, 0)
    r1 = vcat(r1, 0)
    magnitude = norm(cross(r0, r1))
    return magnitude/2
end

"""
get sector area
"""
function getsectorarea(theta1, theta2, rho)
    angle = abs(theta1 - theta2)
    return 1/2*angle*rho^2
end

n = 2
@var x[1:n]
rng = MersenneTwister(2)
nsim = 3
thetavec = collect(0:0.05:2*pi)
for sim in 1:nsim
    for alpha in 0:0.1:1
        #get parameters leading to feasibility
        #r0, A, B = getfeasiblepars(n, alpha, rng)
        r0, A, B = sampleparameters(n, rng)
        pars = (alpha, r0, A, B)
        #solve unperturbed system
        solmat0 = makeandsolve(x, pars)
        #get positions of feasible equilibria
        indrows = getfeasrowinds(solmat0)
        #get the matrix of feasible equilibria
        feas_eq_mat = solmat0[indrows,:]
        #find feasibility domain for each feasible equilibria
        for eq in 1:length(indrows)
            #store equlibrium we are dealing with
            xstari = feas_eq_mat[eq,:]
            #initialize polygon area
            rprev = []
            area = 0
            #traverse full 3D space in polar coordinates 
            rho = 0
            rnew = r0
            feasible = true #start with feasible equilibrium always
            limit = 0 # do I need this?
            while feasible  
                #make a round about
                for i_theta in 1:length(thetavec)
                    theta = thetavec[i_theta]
                    println("Searching boundary for sim = ", sim, " alpha = ", alpha, " rho = ", rho, " theta = ", theta)
                    rnew = r0 .+ pol2cart(theta, rho)
                    #form new parameter set 
                    parsnew = (alpha, rnew, A, B)
                    solmat = makeandsolve(x, parsnew)
                    feasible = anyrowfeasible(solmat)
                    #when we hit a boundary, skip iteration
                    if feasible == false
                        continue
                    end
                    #identify correct equilibria
                    ind_eq_new = identifyequilibrium(xstari, solmat)
                    xstari = solmat[ind_eq_new,:]
                    #check if there is still a feasible equilibria or not
                    feasible = all(xstari .> 0)
                    #add sector area if feasibility remains
                    if i_theta == 1
                        areai = 0
                    else
                        theta1 = thetavec[i_theta-1]
                        theta2 = theta
                        areai = getsectorarea(theta1, theta2, rho)
                    end
                    area += areai
                    storexstar = hcat(sim, alpha, theta, rho, transpose(rnew - r0), limit, transpose(xstari), eq, area)
                    open("../data/boundaryportrait.csv", "a") do io
                        writedlm(io, storexstar, ' ')
                    end
                end
                #increase radius
                rho += 0.1
                #stop simulating when reaching maximum radius
                if rho > 1
                    feasible = false
                    limit = 1
                end
            end
        end
    end
end

# n = 3
# @var x[1:n]

# for alpha in 0:0.5:1
#     #get parameters leading to feasibility
#     r0, A, B = getfeasiblepars(n, alpha, rng)
#     #traverse full 3D space in spherical coordinates 
#     for theta in 0:0.15:2*pi    
#         for phi in 0:0.15:pi
#             rho = 0
#             rnew = r0
#             feasible = true
#             limit = 0
#             #increase magnitud of vector until feasibility is lost
#             while feasible
#                 println("Searching boundary for alpha = ", alpha, "theta = ", theta, " phi = ", phi, " rho = ", rho)
#                 rnew = r0 .+ sph2cart(theta, phi, rho)
#                 #form new parameter set 
#                 parsnew = (alpha, rnew, A, B)
#                 solmat = makeandsolve(x, parsnew)
#                 #check if there is still a feasible equilibria or not
#                 feasible = anyrowfeasible(solmat)
#                 #deal with multiple equilibria in this part, but later###############################
#                 #interrupt search if domain is too big
#                 if rho > 5
#                     feasible = false
#                     limit = 1
#                 end
#                 #increase radius
#                 rho += 0.15 
#             end
#             #record r for which feasibility is lost
#             tosave = hcat(alpha, transpose(rnew-r0), limit)
#             open("../data/feasibility_boundary3spp.csv", "a") do io
#                 writedlm(io, tosave, ' ')
#             end
#         end
#     end
# end
# #plot minimum magnitudes as a function of alpha. Should see an increasing line.


#TODOS

# check why is the basic idea not working?
# go back to the basics of my theory, and try to follow it in the simulations