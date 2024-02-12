using LinearAlgebra
using HomotopyContinuation
using DelimitedFiles
using Random
using IterTools

include("experiment.jl")

function getstorerows(solutionslong, nsols, sim, n, signs, parvalue, eigenvals) 
    nvec = repeat([n], nsols*n) 
    signsmat = repeat(collect(signs)', nsols*n)
    simvec = repeat([sim], nsols*n)
    solcomp = repeat(collect(1:n), nsols)
    eqid = repeat(collect(1:nsols), inner = n)
    eigenvalslong = repeat(eigenvals, inner = n)
    parvalvec = repeat(parvalue, nsols*n)
    return [simvec nvec signsmat parvalvec solcomp eqid solutionslong eigenvalslong]
end

function remove_negative_rows(matrix)
    # Initialize an empty array to store the rows without negative numbers
    filtered_rows = []

    # Iterate through each row in the matrix
    for row in eachrow(matrix)
        # Check if the row contains any negative number
	if all(real(row) .>= 0)
            # If all elements are non-negative, add the row to the filtered_rows array
            push!(filtered_rows, row)
        end
    end

    # Convert the filtered_rows array into a matrix and return it
    return filtered_rows
end

function traversealpha(par_syst, endpars, sim, n, step, signs, x)
    #res = solve(par_syst, initialsol; start_parameters = startpars, 
    #            target_parameters = endpars)	
    #substitute parameters
    syst = parametersubstitution(par_syst, endpars, par_syst.parameters)
    res = solve(syst ./ x, compile = false) #only find solutions with non-zeros
    solmat = real_solutions(res)
    if isempty(sols)
	#no real solutions
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
	tosave = getstorerows(solutionslong, nsols, sim, n, signs, endpars, smallesteigenvals)
	    #save
	    open("/Users/pablolechon/Desktop/pert_hoi/data/portraitsigned.csv", "a") do io
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
	return traversealpha(par_syst, endpars, sim, n, step, signs, x)
    end
end

"""
sample parameters such that ones(n) is a zero of the system
"""
function sampleparameters(n, rng)
    randr = randn(n)
    randA = randn((n,n))
    randB = randn((n,n,n))
    return randr, randA, randB
end

"""
generate all possible sign parametrizations
"""
function allsignedcombs(elements)
    tupvec = product(fill(elements, length(elements))...)
    return tupvec
end

"""
given a number and a parameter, decide if transform it or not
"""
function number2sign(number, parameter)
    if number == 0
	return parameter
    elseif number == 1
	return abs.(parameter)
    else number == -1
	return -abs.(parameter)
    end
end

"""
given a number, sample parameters with appropriate sign parametrization
"""
function samplesignedparametrization(signstuple, n, rng)
    #sample random parameters
    pars = sampleparameters(n, rng)
    signedpars = (number2sign(signstuple[i], pars[i]) for i in 1:3)
    return signedpars
end


"""
main function to run script
"""
function main()
    nmax = 5
    nsim = 200
	    seed = 1 #abs(rand(Int))
    rng = MersenneTwister(seed)
    initialstep = .01
    signedcombs = allsignedcombs([0, 1, -1])
    @var alpha[1:1]
    for n in 2:nmax
	@var x[1:n]
    	for sim in 1:nsim
	    #sample parameters 
	    modelpars = sampleparameters(n, rng)
	    #loop through all possible signed parametrizations
	    for signedcomb in signedcombs
		println("diversity: ", n, " simulation: ", sim, " parametrezation: ", signedcomb)
		#get correct signs on each parameter set
		signedmodelpars = [number2sign(signedcomb[i], modelpars[i]) for i in 1:3]
		pars = (alpha, signedmodelpars...)
		#build system
		eqs = buildglvhoi(pars, x)
		syst = System(eqs, parameters = alpha)
		traversealpha(syst, [.0], sim, n, initialstep, signedcomb, x)
	    end
	end
    end
    return 0
end
main()
