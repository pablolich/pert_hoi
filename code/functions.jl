#############################################################################################
#############################################################################################
#############################################################################################
##########################FUNCTIONS DEALING WITH PARAMETER SAMPLING##########################
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

#############################################################################################
#############################################################################################
#############################################################################################
#########################FUNCTIONS TO DEALING WITH THE GLV HOI MODEL#########################

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

#############################################################################################
#############################################################################################
#############################################################################################
##############################FUNCTIONS TO IDENTIFY CRITICAL RADIUS##########################

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
given a column vector nx1 and a matrix kxn, 
return the row of the matrix that is closest to the vector. 
"""
function closest_row(vector::AbstractVector, matrix::Matrix{Float64})
    # Ensure the vector is a column vector
    if size(vector, 1) != size(matrix, 2)
        throw(ArgumentError("The length of the vector must match the number of columns in the matrix."))
    end

    # Calculate the distances
    distances = [norm(matrix[i, :] .- vector) for i in 1:size(matrix, 1)]

    # Find the index of the minimum distance
    closest_index = argmin(distances)

    # Return the closest row of the matrix
    return matrix[closest_index, :]
end

"""
given all the real solutions to a system of polynomials, and a target vector: 
1. if follow == true, then return the solution that is closest to the target vector
2. if follow == false return all the positive rows.
if the selected row(s) either the closest to target or all positve rows are negative or 
empty, then return nothing, since no positive rows exists.
"""
function select_row(matrix::Matrix{Float64}, 
                    targetsolution::AbstractVector, follow)
    nspp = size(matrix, 2)
    if follow == true
        #select the closest equilibrium to the targetsolution
        return closest_row(targetsolution, matrix)
    else
        selected_rows = []
        # Step 1: Collect all rows with all positive entries
        for row in 1:size(matrix, 1)
            if all(matrix[row, :] .> 0)
                push!(selected_rows, matrix[row, :])
            end
        end
        
        # Step 2: If positive rows exist, return them
        if length(selected_rows) > 0
            return selected_rows #min_elem_row
        else 
            return repeat(-1, nspp)
        end
        
    end
end

"""
perform a disc perturbation of the growthrate in pars of radius rho, and get
the equilibrium responses after perturbation. when interrupt == true the programm will stop
whenever there is a perturbation generating any non-feasible solutions
"""
function perturbondisc(rho, pars, n, nperts, x, interrupt, perturbation_vec, target_equlibria)
    #project a vector of perturbations onto a sphere of radius rho
    perts_rho = rho*perturbation_vec
    equilibria = Matrix{Float64}(undef, 0, n)
    pert = 1
    feasiblesol = true
    while pert < nperts
        #get specific perturbation
        pert_rho_i = perts_rho[pert,:]
        #get the equilibrium achieved for the perturbation with a slightly smaller radius
        target_equilibrium = target_equilibria[pert, :]
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
        if interrupt == true 
            #user would like to interrupt search when it fails to find a feasible equilibrium:
            if nsols == 0 || nfeas_sols == 0 
                #no feasible or real eq for a given perturbation
                println("Condition has been fullfilled with interrupt = ", interrupt, 
                        " and n feas sols = ", nfeas_sols)
                return equilibria
            end
        else
            #user would like to carry all perturbation regardless of sign
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
    end
    return equilibria
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
find middle point between two points
"""
function bisect(x0, x1)
    return 1/2*(x0 + x1)
end

"""
get maximum perturbation on the growth rates retaining feasibility
"""
function findmaxperturbation(rho1, rho2, pars, n, nperts, x, tol,
                             target_equilibria,
                             pertubation_vec, 
                             interrupt=true)
    #initialize number of iterations
    n_iterations = 0
    while abs(rho1 - rho2) > tol
        if n_iterations == 0
            target_equilibria = []
        else 
            target_equilibria = 
        end
        println("[", rho1, ",", rho2, "]")
        #find the middle point of current interval
        global rhob = bisect(rho1, rho2)
        #perform disc perturbation of radius rb
        equilibria = perturbondisc(rhob, pars, n, nperts, x, interrupt)
        #check the minimum x obtained around the disc
        xmin = get_minimum(equilibria, nperts)
        println("xmin = ", xmin)
        #modify interval depending on wehter the minimum is negative or positive
        rho1, rho2 = getlims(rho1, rho2, rhob, xmin)
        #if solution is found, check that no other solutions exist for smaller rs
        if abs(rho1 - rho2) < tol
            for rho in range(rhob, tol, 10)
                equilibria = perturbondisc(rho, pars, n, nperts, x, interrupt)
                xcheck = get_minimum(equilibria, nperts)
                println("Backward checking: ", rho, "x = ", xcheck)
                #recalculate if negative to make sure it's not a mistake of the package (sometimes it happens)
                if xcheck < -tol
                    equilibria = perturbondisc(rho, pars, n, nperts, x, interrupt)
                    xcheck = get_minimum(equilibria, nperts)
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


#############################################################################################
#############################################################################################
#############################################################################################
############################CORRECTED FUNCTIONS FOR CRITICAL RADIUS##########################

"""
a function that takes an nperturbationsxn matrix and determines if all the rows are positive 
(that is, all have only elements greater than 0), and different from each other 
(that is, each row is unique within a certain tolerance). 
If these two conditions are met, return true. otherwise return false.
"""
function check_rows_positive_and_unique(matrix::Matrix{Float64}, tolerance::Float64)
    # Check if all rows are positive
    all_positive = all(all(row .> 0) for row in eachrow(matrix))

    if !all_positive
        return false
    end

    # Check if all rows are unique within the given tolerance
    unique_rows = []

    for row in eachrow(matrix)
        # Check if the current row is similar to any existing unique row
        is_unique = true
        for unique_row in unique_rows
            if all(abs(row[i] - unique_row[i]) <= tolerance for i in 1:length(row))
                is_unique = false
                break
            end
        end

        if is_unique
            push!(unique_rows, row)
        end
    end

    return length(unique_rows) == size(matrix, 1)
end

"""
gets the first set of perturbations. all must be positive, but different (within tolerance)
"""
function get_first_target(parameters::Tuple, 
                          rho::Float64, 
                          perturbations::Matrix{Float64}, 
                          nperturbations::Int64, 
                          x::Vector{Variable},
                          n::Int64)
    all_equilibria = perturbondisc(perturbations, rho, parameters, x, true)
    #confirm that indeed they are different and positive
    first_target = select_equilibria(all_equilibria, 
                                     repeat([1], n),
                                     "positive")
    is_target_valid = check_rows_positive_and_unique(first_target)
    if is_target_valid 
        return first_target
    else
        error("non valid first target. either contains negative elements (tolerance too big)
               or too similar elements (tolerance too small). Check first_target and change
               tolerance accordingly.")
end

"""
once solution is found, run backward checks to confirm it is indeed one
"""
function confirm_solution(rhob::Float64, 
                          pars::Tuple,
                          n::Int64,
                          x::Vector{Variable},
                          nperturbations::Int64,
                          target::Matrix{Float64}, mode::String, 
                          tol::Float64)
    for rho in range(rhob, tol, 10)
        all_equilibria = perturbondisc(rho, pars, n, nperturbations, x, true)
        #some could be non-feasible, if mode == "follow"
        xcheck = get_minimum(all_equilibria, nperturbations, target, mode)
        #println("Backward checking: ", rho, "x = ", xcheck)
        #recalculate if negative to make sure it's not a mistake of the package (sometimes it happens)
        if xcheck < -tol
            equilibria = perturbondisc(rho, pars, n, nperturbations, x, true)
            xcheck = get_minimum(equilibria, nperts)
        end
        #if at some point x becomes negative or complex again, then another 0 exists
        if xcheck < -tol
            return xcheck
    end
    return xcheck
end 

"""
perform a disc perturbation of the growthrate in pars of radius rho, and get
the equilibrium responses after perturbation. when interrupt == true the programm will stop
whenever there is a perturbation generating any non-feasible solutions
"""
function perturbondisc(perturbations::Matrix{Float64},
                       rho::Float64,
                       parameters::Tuple,
                       #n::Int64,
                       x::Vector{Variable},
                       interrupt::Bool)
    equilibria = []
    nperts = size(perturbations, 1)
    #project perturbations to current radius
    perturbations_rho = rho*perturbations
    for i in 1:nperts
        pert_i = pertsurbations_rho[i,:]
        # print(print("\e[A\e[2K"))
        # println("Perturbation number: ", pert,"/", nperts, " n = ", n)
        #form parameters of perturbed system with alpha value
        rpert = parameters[2] + pert_i
        parsnew = (parameters[1], rpert, parameters[3], parameters[4])
        #solve system (get real solutions)
        solmat = makeandsolve(x, parsnew)
        #count solutions (real, feasible)
        nsols = length(solmat)
        nfeas_sols = length(getfeasrowinds(solmat))
        if (nsols == 0 || nfeas_sols == 0) && interrupt
            #end perturbation if desired, whenever feasibility breaks
            return equilibria
        end
        #add equilibria to total
        push!(equilibria, solmat)
    end
    return equilibria
end

"""
Auxiliary function to find the closest row in a matrix to a given vector
"""
function closest_row(vector::Vector{Float64}, matrix::Matrix{Float64})
    min_distance = Inf
    closest_row_index = 0

    # Iterate through each row of the matrix
    for row_index in 1:size(matrix, 1)
        current_row = matrix[row_index, :]
        # Calculate the least squares distance
        distance = sum((current_row .- vector).^2)

        # Update if a closer row is found
        if distance < min_distance
            min_distance = distance
            closest_row_index = row_index
        end
    end

    return closest_row_index
end

"""
Auxiliary function to find a positive row or the closest row if needed
"""
function positive_row(target_vector::Vector{Float64}, matrix::Matrix{Float64})
    positive_rows_indices = []

    # Iterate through each row of the matrix
    for row_index in 1:size(matrix, 1)
        current_row = matrix[row_index, :]
        # Check if all elements in the current row are positive
        if all(current_row .> 0)
            push!(positive_rows_indices, row_index)
        end
    end

    # Select the appropriate row based on the positive rows found
    if length(positive_rows_indices) == 1
        return positive_rows_indices[1]  # Return the only positive row found
    elseif length(positive_rows_indices) > 1
        # If multiple positive rows are found, form a matrix of those rows
        positive_matrix = matrix[positive_rows_indices, :]

        # Use the closest_row function to find the closest row to the target vector
        closest_row_index = closest_row(target_vector, positive_matrix)

        return positive_rows_indices[closest_row_index]  # Return the index from the original matrix
    else
        # If no positive rows found, use the original matrix to find the closest row
        return closest_row(target_vector, matrix)  # Return closest row index from original matrix
    end
end

"""
given a collection of equilibra, and a target, select equilibria based on criteria in mode
if mode == "follow", then the row in each matrix of perturbed equilibria closest to the corresponding
target is selected
if mode == "positive", then a positive equilibrium is selected if the equilibrium resulting form 
follow is non-feasible. 
if, even this mode yields a negative equilibrium, then the followed equilibrium is returned.
"""
function select_equilibria(array_of_matrices::Vector{Matrix{Float64}}, 
                           target_matrix::Matrix{Float64},
                           mode::String)
    k = length(array_of_matrices)   # Number of matrices
    matched_matrix = Matrix{Float64}(undef, k, size(target_matrix, 2))  # Initialize matched matrix
    for i in 1:k
        # Get the ith row of the target matrix
        target_row = target_matrix[i, :]
        if mode == "follow"
            # Find the index of the closest row in the corresponding matrix
            row_index = closest_row(target_row, array_of_matrices[i])
        
        elseif mode == "positive"
            row_index = positive_row(target_row, array_of_matrices[i])
        else
            error("Don't know this mode of select equilibria")
        end
        matched_matrix[i, :] = array_of_matrices[i][row_index, :]
    end
    return matched_matrix
end


"""
given the array of matrices with equilibria from all perturbations, select the minimum of them
after some filtering based on the mode perturbations and result of perturbations
This function is only called inside findmaxperturbation, which assumes interrupt == true always
therefore, the selected_equilibria should always contain positive rows if the mode of selection is
positive. Otherwise the perturbation protocol should have not finished, and all_equilibria should
be shorter than the number of perturbations.
"""
function get_minimum(all_equilibria::Vector{Matrix{T}},
                     nperturbations::Int64,
                     target_equilbria::Matrix{Float64},
                     mode::String
                     )
    if length(all_equilibria) < nperturbations
        #the perturbation protocol stopped early (non-feasibility encountered)
        xmin = -1.0
    else 
        #all perturbations yielded feasible equilibria
        selected_equilibria = select_equilibria(all_equilibria,
                                                target_equilibria,
                                                mode)
        #this minimum could be negative if mode=="follow"
        xmin = minimum(selected_equilibria)
    end
    return xmin
end

"""
get maximum perturbation on the growth rates retaining feasibility
"""
function findmaxperturbation(rho1::Float64, rho2::Float64,
                             perturbations::Matrix{Float64},
                             nperturbations::Int64,
                             parameters::Tuple, 
                             x::Vector{Variable},
                             target_equilibria::Matrix{Float64},
                             mode::String,
                             tol::Float64=1e-9
                             )
    #println("[", rho1, ",", rho2, "]")
    #find the middle point of current interval
    global rhob = bisect(rho1, rho2)
    #perform disc perturbation of radius rb
    while abs(rho1 - rho2) > tol
        all_equilibria = perturbondisc(perturbations, rhob, parameters, x, true)
        #some could be non-feasible, if mode == "follow"
        xmin = get_minimum(all_equilibria, nperturbations, target_equilbria, mode)
        if (xmin > 0)
            #assign new target equilibria if all perturbations are feasible
            target_equilbria = selected_equilibria
        end
        #redefine limits
        rho1, rho2 = getlims(rho1, rho2, rhob, xmin)
        #if solution is found, check that no other solutions exist for smaller rs
        if abs(rho1 - rho2) < tol
            xmin = confirm_solution(rhob, parameters, n, x, nperturbations, target, mode, tol)
            if xmin < -tol
                rho1 = 0 
                rho2 = rhob
            end
        end
        return findmaxperturbation(rho1, rho2, perturbations, nperturbations, parameters,
                                   x, target_equilibria, mode, tol)
    end
    return rhob, target_equilibria
end

"""
function to count the average number of positve rows in each matrix of equilibria arising from 
a perturbation.
"""
function average_number_positive_rows(array_of_matrices::Vector{Matrix{T}})
    k = length(array_of_matrices)   # Number of matrices
    # Initialize vector of number of positive rows in each matrix
    n_positive_rows = []  
    for i in 1:k #loop through each matrix
        matrix_i = array_of_matrices[i]
        # Iterate through each row of the matrix
        for row_index in 1:size(matrix_i, 1)
            current_row = matrix[row_index, :]
            # Check if all elements in the current row are positive
            if all(current_row .> 0)
                push!(positive_rows_indices, row_index)
            end
        end
        push!(n_positive_rows, length(positive_rows_indices))
    end
    return mean(n_positive_rows)
end

"""
function to compute the proportion of feasible states of a given equilibrium matrix
"""
function proportion_of_positive_rows(matrix::Matrix{Float64})
    nrows = size(matrix, 1)
    positive_rows_indices = []
    # Iterate through each row of the matrix
    for row_index in 1:nrows
        current_row = matrix[row_index, :]
        # Check if all elements in the current row are positive
        if all(current_row .> 0)
            push!(positive_rows_indices, row_index)
        end
    end
    return length(positive_rows_indices)/nrows
end