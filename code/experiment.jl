using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Random
using IterTools

"""
evaluate jacobian of glv hoi model and return smallest eigenvalue
"""
function linearstability(syst, point)
    jac = jacobian(syst, point)
    eigRe = real(eigen(jac).values)
    return maximum(eigRe)
end

function linearstabilitymany(syst, points)
    nsols = size(points, 1)
    return [linearstability(syst, points[i]) for i in 1:nsols]
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
perturb growthrates
"""

"""
maximize the timestep that allows you to get a smooth curve on the phase diagram
"""
function getstep(maxdistance, initialstep)
    dist = norm(initialeq - finaleq)
    return parameterrange*(exp(1/dist) - 1)
end


"""
substitute parameter by actual value in system
"""
function parametersubstitution(par_syst, par_value, pars)
    neqs = length(par_syst)
    copysyst = deepcopy(par_syst)
    for i in 1:neqs
 	copysyst.expressions[i] = subs(copysyst.expressions[i] , pars => par_value)
    end
    return copysyst
end

function count_nonzero(matrix::Vector{Vector{ComplexF64}}; tol::Float64=1e-9)
    return [count(x -> abs(real(x)) > tol, row) for row in matrix]
end

"""
auxiliary function to facilitate storing results
"""
function getstorerows(solutionslong, nsols, sim, n, constr, parvalue, eigenvals) 
    nvec = repeat([n], nsols*n) 
    simvec = repeat([sim], nsols*n)
    solcomp = repeat(collect(1:n), nsols)
    eqid = repeat(collect(1:nsols), inner = n)
    eigenvalslong = repeat(eigenvals, inner = n)
    parvalvec = repeat(parvalue, nsols*n)
    constrvec = repeat([constr], nsols*n)
    return [simvec nvec constrvec parvalvec solcomp eqid solutionslong eigenvalslong]
end

"""
decompose solution in real and imaginary, and stack components in columns
"""
function decomposesolution(solution)
    n = length(solution)
    #initialize matrix
    sol_mat = Array{Float64}(undef, n)
    for i in 1:n
        sol_mat[i,1] = real(solution[i])
    end
    return sol_mat
end

"""
perform the previous decomposition for many roots
"""
function decomposemany(solutionsmat)
    nsols = size(solutionsmat, 1)
    nspp = length(solutionsmat[1])
    #initialize matrix of storing many decomposed solutions
    sols_mat = Array{Float64}(undef, nspp*nsols, 1)
    for i in 1:nsols
        solutioni = solutionsmat[i]
        decomposedsol = decomposesolution(solutioni)
        #store
        sols_mat[(nspp*(i-1)+1):(nspp*i),:] = decomposedsol
    end
    return sols_mat
end

"""
sample tesnor of HOIs with constraints
"""
function constrainB(B, r0, n, rng)
    Bconst = zeros((n, n, n))
    for i in 1:n
	randmat = randn(rng, (n,n))
	Bconst[:,:,i] .= -r0/sum(B[:,:,i]) .* B[:,:,i]
    end
    return Bconst
end

"""
sample row stochastic matrix
"""
function constrainA(A, r0)
    sumrows = sum(A, dims = 2)
    Aconst = -r0 .* diagm(1 ./ vec(sumrows)) * A
    return Aconst
end

"""
sample parameters such that ones(n) is a zero of the system
"""
function sampleparameters(n, rng, constrained = false)
    r0 = randn(rng)
    r = repeat([r0], n)
    A = randn((n,n))
    B = randn((n,n,n))
    Aconstr = constrainA(A, r0)
    Bconstr = constrainB(B, r0, n, rng)
    return r, A, B, Aconstr, Bconstr
end

function traversealpha(par_syst, endpars, sim, n, step, x, constr)
    #res = solve(par_syst, initialsol; start_parameters = startpars, 
    #            target_parameters = endpars)	
    #substitute parameters
    syst = parametersubstitution(par_syst, endpars, par_syst.parameters)
    res = solve(syst ./ x, compile = false) #only find solutions with non-zeros
    solmat = real_solutions(res)
    if isempty(solmat)
	#no real soutions
    else
	#solmat = mapreduce(permutedims, vcat, sols)
	#solmat = remove_negative_rows(solmat) #remove nonfeasible sols
	if isempty(solmat)
	    #no positive solutions
	else 
	    nsols = size(solmat, 1)
	    smallesteigenvals = linearstabilitymany(syst, solmat)
	    #get how many species there are in each solution
	    #save result
	    solutionslong = decomposemany(solmat)
	    tosave = getstorerows(solutionslong, nsols, sim, n, constr, endpars, smallesteigenvals)
	    #save
	    open("/Users/pablolechon/Desktop/pert_hoi/data/portrait.csv", "a") do io
		writedlm(io, tosave, ' ')
	    end
	end
    end
    #prepare for next step
    endpars .+= step
    #check if traversing has converged
    if endpars[1] >= 1.0
	return 0
    #run function again if not
    else
	return traversealpha(par_syst, endpars, sim, n, step, x, constr)
    end
end

"""
main function to run script
"""
function main()
    nmax = 5
    nsim = 100
    seed = 1 #abs(rand(Int))
    rng = MersenneTwister(seed)
    initialstep = .01
    @var alpha[1:1]
    for n in 2:nmax
	@var x[1:n]
    	for sim in 1:nsim
	    println("diversity: ", n, " simulation: ", sim)
	    #sample parameters (unconstr)
	    allpars = sampleparameters(n, rng)
	    modelpars = allpars[1:3]
	    modelparsconstr = (allpars[1], allpars[4:5]...)
    	    parsunconstr = (alpha, modelpars...)
	    parsconstr = (alpha, modelparsconstr...)
	    #build system
	    eqsunconstr = buildglvhoi(parsunconstr, x)
	    systunconstr = System(eqsunconstr, parameters = alpha)
	    traversealpha(systunconstr, [.0], sim, n, initialstep, x, false)
	    #recalculate for constrained parameters
	    eqsconstr = buildglvhoi(parsconstr, x)
	    systconstr = System(eqsconstr, parameters = alpha)
	    traversealpha(systconstr, [.0], sim, n, initialstep, x, true)
	end
    end
    return 0
end

main()
