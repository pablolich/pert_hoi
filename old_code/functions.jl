#code to store needed functions for simulations

#load needed packages
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Random
using JuMP
using Ipopt
using Kronecker
using Plots, Colors
using JLD


#functions

"""
sample row stochastic matrix
"""
function constrainA(A, r0)
    sumrows = sum(A, dims = 2)
    Aconst = -r0 .* diagm(1 ./ vec(sumrows)) * A
    return Aconst
end

"""
auxiliary constrain-building functions
"""
function buildrowconstraintmat(n)
    return collect(transpose(kronecker(I(n), repeat([1], n))))
end

function buildcolconstraintmat(n)
    return collect(transpose(kronecker(repeat([1], n), I(n))))
end

"""
get the closest (in the frobenious norm sense) bistochastic matrix to a given one
"""
function getclosestbistochastic(A, n)
    #get matrices of  constraints
    rowconsmat = buildrowconstraintmat(n)
    colconsmat = buildcolconstraintmat(n)
    consvec = repeat([1], n)
    #perform optimization
    model = Model(Ipopt.Optimizer)
    @variable(model, x[1:n^2])
    @objective(model, Min, sqrt(sum((vec(A) .- x).^2)))
    @constraint(model, rowconsmat * x .== consvec)
    @constraint(model, colconsmat * x .== consvec)
    optimize!(model)
    #get result in matrix form
    result = reshape([value(x[i]) for i in 1:n^2], n, n)
    return result
end

"""
sample a tensor where all the slices add up to a constant
"""
function getconstantsumB(B0, n, r0)
    Bconstsum = zeros((n,n,n))
    for i in 1:n
	Bi = B0[:,:, i]
	Bconstsum[:,:,i] .= -r0/n*getclosestbistochastic(Bi, n)
    end
    return Bconstsum
end

"""
sample tesnor of HOIs with constraints
"""
function constrainB(B, r0, n)
    Bconst = zeros((n, n, n))
    for i in 1:n
	Bconst[:,:,i] .= -r0/sum(B[:,:,i]) .* B[:,:,i]
    end
    return Bconst
end

"""
sample parameters with appropriate constraints
"""
function sampleparameters(n, rng, constrain_type)

    r0 = randn(rng)
    r = repeat([r0], n)
    A = randn(rng, (n,n))
    B = randn(rng, (n,n,n))
    Aconstr = constrainA(A, r0) #equilibrium preserving constraints
    if constrain_type == 1
        Bconstr = getconstantsumB(B, n, r0) #equilibrium and tractability constraints.
    else
        Bconstr = constrainB(B, r0, n) #only equilibrium preserving constraints.
    end
    return r, Aconstr, Bconstr
end

"""
generate random points on the surface of a n-dimensional hypersphere of radius rho.
when dimension is 2, the points are evenly distributed
"""
function points_hypersphere(dim::Int, rho::Float64, num_points::Int)
    if dim == 2
        points = zeros(num_points, 2)  # Initialize matrix to store points
        for i in 1:num_points
            theta = 2π * (i - 1) / num_points  # Calculate angle for current point
            points[i, 1] = rho * cos(theta)  # Calculate x coordinate
            points[i, 2] = rho * sin(theta)  # Calculate y coordinate
        end
        return points
    else
        points = randn(num_points, dim)  # Generate random points in dim-dimensional space
        norms = [norm(points[i,:]) for i in 1:num_points]  # Calculate norms of each point
        scaled_points = rho * (points ./ norms)  # Scale points to lie on the surface of the sphere of radius rho
        return scaled_points
    end
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
get matrix with real solutions of system of polynomials
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
get angle between two vectors
"""
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

"""
find middle point between two points
"""
function bisect(x0, x1)
    return 1/2*(x0 + x1)
end

"""
get appropriate limits based on sign of the middle point
"""
function getlims(x0, x1, xb, xmin)
    if xmin > 0
        x0 = xb
    else
        x1 = xb
    end
    return x0, x1
end

"""
given all the real solutions to a system of polynomials, return all the positive rows.
if no positive rows exist
"""
function select_row(matrix::Matrix{Float64}, targetsolution)
    positive_rows = []
    # Step 1: Collect all rows with all positive entries
    for row in 1:size(matrix, 1)
        if all(matrix[row, :] .> 0)
            push!(positive_rows, matrix[row, :])
        end
    end
    
    # Step 2: If positive rows exist, return them
    if length(positive_rows) > 0
        return positive_rows #min_elem_row
    end
    
    # Step 3: Compute the norm 2 distance to the vector of ones if no all-positive row found
    ones_vector = ones(size(matrix, 2))
    norms = [norm(matrix[row, :] - ones_vector) for row in 1:size(matrix, 1)]
    
    # Step 4: Find the row with the minimum norm 2 distance
    min_index = argmin(norms)
    
    return [matrix[min_index, :]]
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
perform a disc perturbation of the growthrate in pars of radius rho, and get
the equilibrium responses after perturbation
"""
function perturbondisc(rho, pars, n, nperts, x, interrupt)
    #generate all perturbations on surface of hypershpere of radius rho
    perts_rho = points_hypersphere(n, rho, nperts)
    equilibria = Matrix{Float64}(undef, 0, n)
    pert = 1
    feasiblesol = true
    while pert < nperts
        #get specific perturbation
        pert_rho_i = perts_rho[pert,:]
        #increase counter
        pert += 1
        print(print("\e[A\e[2K"))
        println("Perturbation number: ", pert,"/", nperts, " n = ", n)
        #form parameters of perturbed system with alpha value
        rpert = pars[2] + pert_rho_i
        parsnew = (pars[1], rpert, pars[3], pars[4])
        #solve system (get real solutions)
        solmat = makeandsolve(x, parsnew)
        nsols = length(solmat)
        nfeas_sols = length(getfeasrowinds(solmat))
        if nsols == 0 #only complex solutions
            equilibria = vcat(equilibria, repeat([-Inf], n)') #treat complex solutions as negative.
            feasiblesol = false #we look for radii subtending only positive equilibria
        elseif nfeas_sols == 0 #no feasible solutions
            equilibrium = select_row(solmat)
            equilibria = vcat(equilibria, equilibrium[1]')
            feasiblesol = false
        else #one or more feasible solutions
            #get xmin from selected equilibrium
            equilibrium = select_row(solmat)
            nsols = length(equilibrium)
            for i in 1:nsols
                equilibriumi = equilibrium[i]
                equilibria = vcat(equilibria, equilibriumi')
            end
        end
        if interrupt && feasiblesol == false #no feasible eq for a given perturbation
            println("Condition has been fullfilled with interrupt = ", interrupt, " and feasible sol = ", feasiblesol)
            return equilibria
        end
    end
    return equilibria
end

"""
"""
function plotperturbations(matrix)
    #transform infinties to 0s
    # Get the number of rows in the matrix
    num_rows = size(matrix, 1)
    # Iterate over each row
    for i in 1:num_rows
        # Check if all elements in the row are -Inf
        if all(x -> x == -Inf, matrix[i, :])
            # Replace the row with zeros
            matrix[i, :] .= 0
        end
    end
    scatter(matrix[:,1], matrix[:, 2])
end

"""
get maximum perturbation on the growth rates retaining feasibility
"""
function findmaxperturbation(rho1, rho2, pars, n, nperts, x, tol, interrupt=true)
    while abs(rho1 - rho2) > tol
        println("[", rho1, ",", rho2, "]")
        #find the middle point of current interval
        global rhob = bisect(rho1, rho2)
        #perform disc perturbation of radius rb, and check the minimum x obtained around the disc
        xmin = minimum(perturbondisc(rhob, pars, n, nperts, x, interrupt))
        println("xmin = ", xmin)
        #modify interval depending on wehter the minimum is negative or positive
        rho1, rho2 = getlims(rho1, rho2, rhob, xmin)
        #if solution is found, check that no other solutions exist for smaller rs
        if abs(rho1 - rho2) < tol
            for rho in range(rhob, tol, 10)
                xcheck = minimum(perturbondisc(rho, pars, n, nperts, x, interrupt))
                println("Backward checking: ", rho, "x = ", xcheck)
                #recalculate if negative to make sure it's not a mistake of the package (sometimes it happens)
                if xcheck < -tol
                    xcheck = minimum(perturbondisc(rho, pars, n, nperts, x, interrupt))
                end
                #if at some point x becomes negative or complex again, then another 0 exists
                if xcheck < -tol || xcheck == -Inf
                    rho1 = 0
                    rho2 = rho  
                    break
                end
            end
        end
        return findmaxperturbation(rho1, rho2, pars, n, nperts, x, tol)
    end
    return rhob
end

# npertbase = 10
# n = 2
# nperts = npertbase^n
# #define variables for polynomial construction
# @var x[1:n]
# rng = MersenneTwister(2)
# constrain_type = 1
# # for i in 1:42
# #     sampleparameters(n, rng, constrain_type)
# # end
# r0, A, B = sampleparameters(n, rng, constrain_type)
# alpha = 0.1
# pars = (alpha, r0, A, B)
# # #parsloaded = load("../data/sim43pars.jld")
# # #pars = (alpha, parsloaded["r0"], parsloaded["A"], parsloaded["B"])
# rmax = findmaxperturbation(0, 30, pars, n, nperts, x, 1e-9)
# #rmax = 16
# #checkresultvisually(rmax, pars, n, nperts, x)
# plotperturbations(perturbondisc(rmax, pars, n, nperts, x))

"""
function to analyze the proportion of feasible states as the radius increases from 1 to rinf
"""
function proportion_of_positive_rows(matrix::AbstractMatrix{T}) where T
    # Ensure the matrix is not empty
    if size(matrix, 1) == 0
        return 0.0
    end
    
    # Calculate the number of rows where all elements are positive
    num_positive_rows = length(getfeasrowinds(matrix))
    
    # Calculate the proportion
    proportion = num_positive_rows / size(matrix, 1)
    
    return proportion
end
