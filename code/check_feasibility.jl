#write code to see whats the maximum fully feasible region

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
        A = randn(rng, (n,n))#/n
        B = randn(rng, (n,n,n))#/n^2
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


function angle(v1::Vector{T}, v2::Vector{T}) where T
    dot_product = dot(v1, v2)
    norm_v1 = norm(v1)
    norm_v2 = norm(v2)
    
    # Ensure vectors are not zero vectors
    if norm_v1 == 0.0 || norm_v2 == 0.0
        throw(ArgumentError("Input vectors cannot be zero vectors"))
    end
    
    cos_angle = dot_product / (norm_v1 * norm_v2)
    
    # Ensure cos_angle is within [-1, 1]
    cos_angle = min(max(cos_angle, -1.0), 1.0)
    
    angle_rad = acos(cos_angle)
    return angle_rad
end

function linearresponse(x, pars, x0, r0)
    #make and solve new system
    syst = System(buildglvhoi(pars, x) ./ x)
    #calculate linear approximation of perturbation
    jacmat = jacobian(syst, x0)
    inv_j = []
    try
        inv_j = inv(jacmat)
    catch e
        if isa(e, SingularException)
            println("Matrix is not invertible, continuing...")
            inv_j = zeros(2, 2)
        end
    end
    deltar = pars[2] - r0
    deltax = -inv_j*deltar
    return deltax + x0
end

function generate_limits(dim::Int, lower::Real=0.0, upper::Real=1.0)
    limits = [(lower, upper) for _ in 1:dim]
    return limits
end

function generate_grid(limits::Vector{Tuple{T,T}}, grain::Int64=3) where T
    ranges = [range(limit[1], limit[2], length=grain) for limit in limits]
    points = vcat(collect(Iterators.product(ranges...))...)
    npoints = size(points)[1]
    return [collect(points[i]) for i in 1:npoints]
end

function points_hypersphere(dim::Int, rho::Float64, num_points::Int)
    points = randn(num_points, dim)  # Generate random points in dim-dimensional space
    norms = [norm(points[i,:]) for i in 1:num_points]  # Calculate norms of each point
    scaled_points = rho * (points ./ norms)  # Scale points to lie on the surface of the sphere of radius rho
    return scaled_points
end


n = 2 #start with n = 2 to recover previous results
@var x[1:n]
rng = MersenneTwister(2)
nsim = 43
idrows = 0
#generate all perturbations 
#exponentially distributed alphas
alphavec = exp.(collect(-3:0.3:0))
rhovec = 10 ./ collect(1:1:1000)
nperts = 200
#loop over simulations
for sim in 1:nsim
   #sample parameters
   #get parameters leading to feasibility
   #constr = false
   #alpha = 0  
   #r0, A, B = getfeasiblepars(n, alpha, rng)
   constr = true
   r0, A, B = sampleparameters(n, rng)
   if sim != 43
    continue
   end
   #loop over hypersphere radii
   for alpha in collect(0:0.1:1)
        for rho in rhovec
            #generate all perturbations on surface of hypershp of radius rho
            rperts = points_hypersphere(n, rho, nperts)
            #for each alpha, traverse all perturbations on the hypershpere
            full_feasible = false
            println("Searching boundary for simulation = ", sim, " alhpa = ", alpha, " rho: ", rho)
            for pert in 1:nperts
                #build new r vector
                perti = rperts[pert,:]
                rnew = r0 + perti
                global idrows += 1
                parsnew = (alpha, rnew, A, B)
                #solve perturbed system
                solmat = makeandsolve(x, parsnew)
                nsols = size(solmat)[1]
                if size(solmat)[1] == 0 #no real solutions
                    println("not real solutions, breaking")
                    break
                elseif !anyrowfeasible(solmat) #no feasible solutions
                    println("not feasible solutions, breaking")
                    break
                else
                    simvec = repeat([sim], n*nsols)
                    nvec = repeat([n], n*nsols)
                    alphavec = repeat([alpha], n*nsols)
                    idrowvec = repeat([idrows], n*nsols)
                    sppid = repeat(collect(1:n), nsols)
                    eqid = repeat(collect(1:nsols), inner = n)
                    rvecs = repeat(rnew - r0, nsols)
                    eqvec = vcat(transpose(solmat)...)
                    npertvec = repeat([pert], nsols*n)
                    rhovec = repeat([rho], nsols*n)
                    tosave = hcat(simvec, nvec, alphavec, idrowvec, rhovec, eqid, rvecs, sppid, eqvec, npertvec)
                    open("../data/maximum_fully_feasible43.csv", "a") do io
                        writedlm(io, tosave, ' ')
                    end
                end
                if pert == nperts
                    println("found fully feasible, breaking")
                    full_feasible = true
                end
            end
            if full_feasible == true
                break
            end
       end
   end
end