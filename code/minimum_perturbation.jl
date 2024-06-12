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
	Bconstr = constrainB(B, r0, n) #only equilibrium preserving constraints.
    #Bconstrstoch = getconstantsumB(B, n, r0) #equilibrium and tractability constraints.
    return r, Aconstr, Bconstr
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

function getpermutations(n)
    nperm = 2^n-1
    mat = zeros(nperm, n)
    for i in 1:nperm
        perm_i = digits(i, base = 2, pad = n)
        mat[i,:] = perm_i
    end
    return mat
end

function cartesianperturbations(vec, allperts, n)
    #get all permutations of n elements
    nperm = 2^n-1
    nperts = length(allperts)
    perm_mat = getpermutations(n)
    pert_mat = []
    #for each permutation (row of perm_mat) form all perturbations
    for i in 1:nperm
        pert_mask = perm_mat[i,:]
        pert_mat_i = [ vec .+ allperts[i] .* pert_mask for i in 1:nperts ]
        println(pert_mat_i)
        push!(pert_mat, pert_mat_i) 
    end
    return vcat(pert_mat...)
end

function generate_limits(dim::Int, lower::Real=0.0, upper::Real=1.0)
    limits = [(lower, upper) for _ in 1:dim]
    return limits
end

function generate_grid(limits::Vector{Tuple{T,T}}, grain::Int64=3) where T
    ranges = [range(limit[1], limit[2], length=grain) for limit in limits]
    points = collect(Iterators.product(ranges...))
    return points
end

n = 2
@var x[1:n]
rng = MersenneTwister(2)
nsim = 200
thetavec = collect(0:pi/64:2*pi)
for sim in 1:nsim
    #get parameters leading to feasibility
    constr = false
    alpha = 0
    r0, A, B = getfeasiblepars(n, alpha, rng)
    #constr = true
    #r0, A, B = sampleparameters(n, rng)
    for alpha in collect(0:0.1:1)
        pars = (alpha, r0, A, B)
        #solve unperturbed system
        solmat0 = makeandsolve(x, pars)
        #get positions of feasible equilibria
        indrows = getfeasrowinds(solmat0)
        #get the matrix of feasible equilibria
        feas_eq_mat = solmat0[indrows,:]
        if size(feas_eq_mat)[1] == 0
            continue
        end
        #pick reference equilibrium
        if constr
            xstar0 = [1,1]
        else
            xstar0 = feas_eq_mat[1,:]
        end
        #store equlibrium we are dealing with
        xstari = [1,1]
        #initialize polygon area
        rprev = []
        area = 0
        #traverse full 3D space in polar coordinates 
        rho = 0
        rnew = r0
        feasible = true #start with feasible equilibrium always
        while feasible  
            #make a round about
            for i_theta in 1:length(thetavec)
                theta = thetavec[i_theta]
                println("Searching boundary for sim = ", sim, " alpha = ", alpha, " rho = ", rho, " theta = ", theta)
                rnew = r0 .+ pol2cart(theta, rho)
                #form new parameter set 
                parsnew = (alpha, rnew, A, B)
                #get linear approximation of response to perturbation
                xstarlinear = linearresponse(x, parsnew, xstar0, r0)
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
                println(angle(rnew, r0))
                storexstar = hcat(sim, alpha, abs(pi/4 - theta),#angle(r0, pol2cart(theta, 1)), 
                                    theta, 
                                    rho, transpose(abs.(r0)), transpose(rnew - r0), transpose(xstari), 
                                    transpose(xstar0), transpose(xstarlinear), area)
                open("../data/boundaryportrait.csv", "a") do io
                    writedlm(io, storexstar, ' ')
                end
            end
            #increase radius
            rho += 0.1
            #stop simulating when reaching maximum radius
            if rho > 1
                feasible = false
            end
        end
    end
end


n = 3
@var x[1:n]
nsim = 2

for sim in 1:nsim
    #sample parameters
    #get parameters leading to feasibility
    constr = false
    alpha = 0
    r0, A, B = getfeasiblepars(n, alpha, rng)
    #constr = true
    #r0, A, B = sampleparameters(n, rng)
    for alpha in 0:0.5:1
        pars = (alpha, r0, A, B)
        #solve unperturbed system
        solmat0 = makeandsolve(x, pars)
        #get positions of feasible equilibria
        indrows = getfeasrowinds(solmat0)
        #get the matrix of feasible equilibria
        feas_eq_mat = solmat0[indrows,:]
        if size(feas_eq_mat)[1] == 0
            continue
        end
        #pick reference equilibrium
        if constr
            xstar0 = [1,1]
        else
            xstar0 = feas_eq_mat[1,:]
        end
        #store equlibrium we are dealing with
        xstari = repeat([1], n)
        #initialize polygon area
        area = 0
        #traverse full 3D space in polar coordinates 
        rho = 0
        rnew = r0
        feasible = true #start with feasible equilibrium always
        while feasible  
            #make a round about
            for i_theta in 1:length(thetavec)

#         for phi in 0:0.15:pi
#             rho = 0
#             rnew = r0
#             feasible = true
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
#                 end
#                 #increase radius
#                 rho += 0.15 
#             end
#             #record r for which feasibility is lost
#             tosave = hcat(alpha, transpose(rnew-r0))
#             open("../data/feasibility_boundary3spp.csv", "a") do io
#                 writedlm(io, tosave, ' ')
#             end
#         end
#     end
    end
end
# #plot minimum magnitudes as a function of alpha. Should see an increasing line.


#TODOS

####################################################################################

#sometimes I get unfeasible equilibria in my phase portraits, would be nice to remove


####################################################################################

# once recovered, do perturbations that are almost homogeneous, and see
    #specifically, do symmetric perturbations around a small cone
    #since I need to know what is the new equilibrium, I need to follow the thread.
    #Where does the new equilibrium goes after perturbation of r, record that.
# how robustness decreases (or increases...)

    #note that in this section I should change the part of the code that
    #identifies the equilibria with the same components, since the components
    #will not be the same any more. however, they will be almost the same
    #thus, I should identify the equilibria with components most similar to one
    #another

# once I have an idea of what the perturbations do, and how much they can
# depart polarly, then go back to the real shit, and draw the allowed 
# angle on r, and see how it translates to the equilibria.\
# if I color the points accordingly to "almost homogeneous perturbations, 
# then I should see that as I increase alhpha, the x's become more robust.\
# that would be the result of the paper, and then I can start writing it.

#also, in the paper I should mention that homogeneous perturbations indeed don't do
#anything in the case of linear perturbations, but this is not the case when 
#there exist Higher-Order 


######################################################################################

#shit works

#to confirm: get the distances to each point, and associated with delta theta. show
# that as alpha increases, the delta theta becomes smaller (that is)

#nice plot search. Plot panes for all simulations in a pdf and pick an nice one to show 
# as an example in the paper.

#calculate the average distances of xmags for each rho, as a function of alpha, and plot. 

#23 41 50 53 102 104! 179 203

######################################################################################

#I made my point, I can write a paper.

#my results at least also hold without the tractability constraints
#they don't hold without any constrains, but thats may be because the equilibrium also
#changes.
#make the equilibrium be 1, and see my results

######################################################################################


#I have fixed the bug! now things make a lot of sense.

#come up with a way to fix the fact that there exist no real equilibria.

#Either only analyze parametrizations that always have feasible equilibria, or 
#accept that alpha might destroy feasibility.

#####################################################################################

#why is my plot of alpha versus avmag not matching with linear approximation for low r?


#####################################################################################


#TODOS
#   2. Write code to perform the same task with more than 2 species. 
#   3. Then drop assumptions sequentially for these cases. Do 5 and 10 species.
#   4. Do the whole simulations again, but this time the exploration only happens
#      within the cone. Whenever I hit unfeasible, or established boundary, then stop\
#      and record the max r for which feasibility is lost.
#   5. This is scalable to more than 2 spp, since I just have to do |r| in n dimensions.
#   6. Average over simulations