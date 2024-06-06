using HomotopyContinuation
using LinearAlgebra
using Random
using Plots

include("experiment.jl")

"""
calculate inverse of jacobian
"""
function jac(n, alpha, A, B, point)
    jac = zeros((n,n))
    for i in 1:n
    	for j in 1:n
	    jac[i, j] = (1-alpha) * A[i, j] + alpha * sum((B[:, j, i] + B[j,:, i]).*point)
	end
    end
    return jac
end

"""
sample parameters such that ones(n) is a zero of the system
"""
function sampleparameters(n, rng)
    r0 = randn(rng)
    r = repeat([r0], n)
    A = randn((n,n))
    B = randn((n,n,n))
    Aconstr = constrainA(A, r0) #equilibrium prererving constraints
	Bconstr = constrainB(B, r0, n) #only equilibrium preserving constraints.
    #Bconstrstoch = getconstantsumB(B, n, r0) #equilibrium and tractability constraints.
    return r, Aconstr, Bconstr
end

"""
eigenvalue formula
"""
function getfirstsingval(r0, n, alpha)
    return 1/(r0*(1+alpha))
end

function find_same_column(matrix)
    num_cols = size(matrix, 2)
    matrix = round.(matrix, digits=5)
    for col in 1:num_cols
        # Check if all elements in the column are the same
        if all(x -> x == matrix[1, col], matrix[:, col])
            return col
        end
    end
    # If no such column found, return -1 or any other suitable indicator
    return -1
end

function findnegativecolumns(matrix::Matrix{Float64})
    num_cols = size(matrix, 2)
    negative_cols = Int[]

    for j in 1:num_cols
        if any(x -> x < 0, matrix[:, j])
            push!(negative_cols, j)
        end
    end

    return negative_cols
end

function find_column(vector::Vector{Float64}, matrix::Matrix{Float64})
    n, m = size(matrix)
    for j in 1:m
	if norm(vector .- matrix[:, j], 2) < 1e-6
            return j
        end
    end
    return false  # Vector not found in any column
end


function identifyequilibria(v::Vector{Vector{Float64}}, tol::Float64=1e-6)
    # Transform the vector of vectors into a matrix
    matori = hcat(v...)
    mat = hcat(v...)

    # Find columns with at least one element less than zero
    cols_to_delete = findnegativecolumns(mat)
    if cols_to_delete != false
	# Delete columns with at least one element less than zero
	matdel = mat[:, setdiff(1:size(mat, 2), cols_to_delete)]
    else
	matdel = mat
    end

    # Find columns where all elements are the same up to a tolerance
    for i in 1:size(matdel, 2)
        col = matdel[:, i]
        if all(x -> abs(x - col[1]) ≤ tol, col)
	    #return col
	    return find_column(col, matori)
        end
    end

    return false  # If no such column is found
end

# nsims = 2000
# nmax = 4
# @var alpha[1:1]
# seed = 1 #abs(rand(Int))
# rng = MersenneTwister(seed)
# pertvals = [0.01, 0.1, 1]
# for n in 2:nmax
#     @var x[1:n]
#     #sample parameters and perturb rs
#     modelpars = sampleparameters(n, rng)
#     r0 = modelpars[1][1]
#     allpars = (alpha, modelpars...)
#     #set up models
#     eqs = buildglvhoi(allpars, x)
#     par_syst = System(eqs ./ x, parameters = alpha)
#     #traverse through alphas
#     for alphavalue in 0:0.1:1
# 		alphavalue = [alphavalue]
# 		#get system at alpha
# 		syst = parametersubstitution(par_syst, alphavalue, par_syst.parameters)
# 		#HERE I CAN IMPLEMENT INTERRUPTION WHEN FINDING THE SOLUTION WITH ALL ONES
# 		#solve system and get real solutions
# 		res = solve(syst, compile = false)
# 		solmat = real_solutions(res)
# 		#for each alpha, do a bunch of perturbations
# 		for pertmod in pertvals
# 			for sim in 1:nsims
# 				println("n: ", n, " alpha: ", alphavalue[1], " pert magnitude: ", pertmod,  " simulation: ", sim)
# 				#create and solve a perturbation of the model
# 				modelparspert = deepcopy(modelpars)
# 				modelparspert[1] .= modelparspert[1]*(1+pertmod*randn())
# 				rpert = modelparspert[1][1]
# 				#calculate linear approximation of perturbation
# 				#originalsyst = System(diagm(x) * syst.expressions)
# 				jacmat = jacobian(syst, ones(n))
# 				inv_j = []
# 				try
# 					inv_j = inv(jacmat)
# 				catch e
# 					if isa(e, SingularException)
# 					println("Matrix is not invertible, continuing...")
# 					continue
# 					else
# 					end
# 				end
# 				deltar = modelparspert[1] - modelpars[1]
# 				deltax = -inv_j*deltar
# 				#create model
# 				allparspert = (alpha, modelparspert...)
# 				eqspert = buildglvhoi(allparspert, x)
# 				par_systpert = System(eqspert ./ x, parameters = alpha)
# 				systpert = parametersubstitution(par_systpert, alphavalue, par_systpert.parameters)
# 				#solve perturbation and get real solutions only
# 				respert = solve(systpert, compile = false)
# 				nrealsols = nreal(respert)
# 				if nrealsols == 0
# 					continue
# 				end
# 				solmatpert = real_solutions(respert)
# 				#compare equilibria
# 				indeq = identifyequilibria(solmatpert) #find equilibria of all species the same  
# 				if indeq != false
# 					xpert = solmatpert[indeq]
# 					deltaxnorm = norm(xpert - ones(n), 2)
# 					else
# 					deltaxnorm = -1
# 				end
# 				#return allpars, allparspert, syst, systpert, solmatpert, xpert, deltax, deltar
# 				tostore = [n alphavalue pertmod sim norm(deltar, 2) norm(deltax, 2) deltaxnorm]
# 				#save
# 				open("/Users/pablolechon/Desktop/pert_hoi/data/perturbationdroptractpres.csv", "a") do io
# 					writedlm(io, tostore, ' ')
# 				end
# 			end
# 		end
#     end
# end
