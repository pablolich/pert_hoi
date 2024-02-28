using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Random
using JuMP
using Ipopt
using Kronecker

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
function constrainB(B, r0, n)
    Bconst = zeros((n, n, n))
    for i in 1:n
	Bconst[:,:,i] .= -r0/sum(B[:,:,i]) .* B[:,:,i]
    end
    return Bconst
end

function isstoc(vec, tol)
    return all(abs.(vec .- 1) .< tol)
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
function sampleparameters(n, rng)
    r0 = randn(rng)
    r = repeat([r0], n)
    A = randn((n,n))
    B = randn((n,n,n))
    Aconstr = constrainA(A, r0)
    Bconstr = constrainB(B, r0, n)
    return r, A, B, Aconstr, Bconstr
end

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

function getconstantsumB(B0, n, r0)
    Bconstsum = zeros((n,n,n))
    for i in 1:n
	Bi = B0[:,:, i]
	Bconstsum[:,:,i] .= -r0/n*getclosestbistochastic(Bi, n)
    end
    return Bconstsum
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
    nmax = 4
    nsim = 500
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
	    parsconstr2 = deepcopy(parsconstr)
	    #get constant sum Bo
	    println(size(allpars[3]))
	    Bconst = getconstantsumB(allpars[3], n, allpars[1][1]) 
	    parsconstr2[4] .= Bconst
	    #build system
	    eqsunconstr = buildglvhoi(parsunconstr, x)
	    systunconstr = System(eqsunconstr, parameters = alpha)
	    traversealpha(systunconstr, [.0], sim, n, initialstep, x, 0)
	    #recalculate for parameters with planted eq
	    eqsconstr = buildglvhoi(parsconstr, x)
	    systconstr = System(eqsconstr, parameters = alpha)
	    traversealpha(systconstr, [.0], sim, n, initialstep, x, 1)
	    #recalculate for constrained parameters  
	    eqsconstr2 = buildglvhoi(parsconstr2, x)
	    systconstr2 = System(eqsconstr2, parameters = alpha)
	    traversealpha(systconstr2, [.0], sim, n, initialstep, x, 2)
	end
    end
    return 0
end

#main()
