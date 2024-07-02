#code to store needed functions for simulations

#load needed packages
using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Random
using JuMP
using Ipopt
using Kronecker


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
generate random points on the surface of a n-dimensional hypersphere of radius rho
"""
function points_hypersphere(dim::Int, rho::Float64, num_points::Int)
    points = randn(num_points, dim)  # Generate random points in dim-dimensional space
    norms = [norm(points[i,:]) for i in 1:num_points]  # Calculate norms of each point
    scaled_points = rho * (points ./ norms)  # Scale points to lie on the surface of the sphere of radius rho
    return scaled_points
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