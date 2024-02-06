using HomotopyContinuation
using RandomMatrix
using LinearAlgebra
using DelimitedFiles
using Random

"""
build glv polynomial
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
    return eqs
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

"""
calculate equilibrium of glv as a function of parameters (r, A, B)
and using parameter homotopy for alpha, which takes values tarpars to endpars

par_syst:
initialsol:
startpars:
endpars:
"""
function traversealpha(par_syst, initialsol, startpars, endpars, sim, n, step)
    #res = solve(par_syst, initialsol; start_parameters = startpars, 
    #            target_parameters = endpars)	
    #substitute parameters
    syst = parametersubstitution(par_syst, endpars, par_syst.parameters)
    res = solve(syst, compile = false)
    solmat = solutions(res)
    nsols = nsolutions(res)
    #save result
    solutionslong = decomposemany(solmat)
    tosave = getstorerows(solutionslong, nsols, sim, n, endpars)
    #save
    open("/Users/pablolechon/Desktop/hoi_pert/data/solutionevoloution.csv", "a") do io
	writedlm(io, tosave, ' ')
    end
    #prepare for next step
    initialsol = solmat[1]
    startpars .= endpars
    endpars .+= step
    #check if traversing has converged
    if startpars[1] >= 1.0
	return 0
    #run function again if not
    else
	return traversealpha(par_syst, initialsol, startpars, endpars, sim, n, step)
    end
end

"""
auxiliary function to facilitate storing results
"""
function getstorerows(solutionslong, nsols, sim, n, parvalue) 
    nvec = repeat([n], nsols*n) 
    simvec = repeat([sim], nsols*n)
    solcomp = repeat(collect(1:n), nsols)
    parvalvec = repeat(parvalue, nsols*n)
    return [simvec nvec parvalvec solcomp solutionslong]
end

"""
decompose solution in real and imaginary, and stack components in columns
"""
function decomposesolution(solution)
    n = length(solution)
    #initialize matrix
    sol_mat = Array{Float64}(undef, n,2)
    for i in 1:n
        sol_mat[i,1] = real(solution[i])
        sol_mat[i,2] = imag(solution[i])
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
    sols_mat = Array{Float64}(undef, nspp*nsols, 2)
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
function sampleB(n, r0, rng)
    B = zeros((n, n, n))
    for i in 1:n
	randmat = randn(rng, (n,n))
	B[:,:,i] .= -r0/sum(randmat) .* randmat
    end
    return B
end

"""
sample row stochastic matrix
"""
function sampleA(n, r0, rng)
    randmat = randn(rng, (n,n))
    sumrows = sum(randmat, dims = 2)
    A = -r0 .* diagm(1 ./ vec(sumrows)) * randmat
    return A
end

"""
sample parameters such that ones(n) is a zero of the system
"""
function sampleparameters(n, rng)
    r0 = randn(rng)
    r = repeat([r0], n)
    A = sampleA(n, r0, rng)
    B = sampleB(n, r0, rng)
    return r, A, B
end

"""
main function to run script
"""
function main()
    nmax = 2
    nsim = 1
    seed = 2
    rng = MersenneTwister(seed)
    initialstep = .01
    @var alpha[1:1]
    for n in 2:nmax
	@var x[1:n]
    	for sim in 1:nsim
	    println("diversity: ", n, " simulation: ", sim)
	    #sample parameters 
	    modelpars = sampleparameters(n, rng)
    	    pars = (alpha, modelpars...)
	    #build system
	    eqs = buildglvhoi(pars, x)
	    syst = System(eqs, parameters = alpha)
	    traversealpha(syst, ones(n), [.0], [.0+initialstep], sim, n, initialstep)
	end
    end
    return 0
end

main()
