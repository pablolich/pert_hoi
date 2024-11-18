using HomotopyContinuation
using LinearAlgebra
using DelimitedFiles
using Statistics #calculate mean()
using Random
using JuMP
using Ipopt
using Kronecker
using Plots, Colors

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
    # Get matrices of constraints
    rowconsmat = buildrowconstraintmat(n)
    colconsmat = buildcolconstraintmat(n)
    consvec = repeat([1.0], n)  # Ensure it's a Float vector for consistency

    # Perform optimization with suppressed output
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "print_level", 0)  # Suppress output
    @variable(model, x[1:n^2])
    @objective(model, Min, sqrt(sum((vec(A) .- x).^2)))
    @constraint(model, rowconsmat * x .== consvec)
    @constraint(model, colconsmat * x .== consvec)

    optimize!(model)

    # Get result in matrix form
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
function buildglvhoi(pars, x, symb_pars::Bool=false)
    n = length(pars[2])
    if symb_pars == false
        #unpack parameters
        alpha, r, A, B = pars
    else
        alpha, _, A, B = pars
        @var r[1:n]
    end
    eqs = r + (1 .- alpha) .* A*x
    #add HOIs
    for i in 1:n
        eqs[i] += (alpha .* ( x'*B[:,:,i]*x ))[1]
    end
    return diagm(x) * eqs
end

function has_row_all_ones(matrix::Matrix{Float64}, tol::Float64=1e-6)
    # Check each row of the matrix
    for row in eachrow(matrix)
        if all(abs.(row .- 1.0) .< tol)
            return true  # Found a row of all ones
        end
    end
    return false  # No row of all ones found
end

"""
get matrix with real solutions of system of polynomials
"""
function makeandsolve(x::Vector{Variable}, 
                      pars::Tuple, 
                      n::Int64, 
                      perturbation::Vector{Float64},
                      only_real_sols::Bool=true,
                      mode::String="all")
    if mode == "all"
        # Form parameters of perturbed system with alpha value
        rpert = pars[2] + perturbation
        parsnew = (pars[1], rpert, pars[3], pars[4])
        #make and solve new system, getting all solutions (or all real solutions)
        syst = System(buildglvhoi(parsnew, x))
        res = solve(syst ./ x, compile = false)
    elseif mode == "follow"
        #solve the system by using parameter homotopy using the target equilibrium as initial point
        @var r[1:n]
        syst = System(buildglvhoi(pars, x, true), parameters = r)
        start_solutions = [ones(n)]
        res = solve(syst, start_solutions; 
                    start_parameters = pars[2],
                    target_parameters = pars[2] .+ perturbation)
    end
    #pick which solutions to return
    if only_real_sols
        solvecs = real_solutions(res)
        if nreal(res) == 0
            solmat = Matrix{Float64}(undef, 0, n)
        else
            solmat = mapreduce(permutedims, vcat, solvecs)
        end
    else
        #return all roots
        solvecs = solutions(res)
        if isempty(solvecs) #solution became complex, but parameter homotopy could not follow it
            solmat = Matrix{ComplexF64}(undef, 0, n)
        else
            solmat = mapreduce(permutedims, vcat, solvecs)
        end
    end
   return solmat
end

"""
Get indices where there is a feasible equilibrium.
The function can handle matrices of real or complex numbers.
"""
function getfeasrowinds(matrix::Matrix{T}, tol::Real = 1e-10) where T
    rowinds = Int[]  # Initialize an empty array for indices
    rowind = 1
    
    for row in eachrow(matrix)
        # Check if the first element is complex
        if typeof(row[1]) <: Complex
            # If any imaginary part is above the tolerance, skip this row
            if any(imag(x) > tol for x in row)
                rowind += 1
                continue
            end
            # If all imaginary parts are below tolerance, check real parts
            if all(real(x) > 0 for x in row)
                push!(rowinds, rowind)
            end
        else
            # For real numbers, proceed as before
            if all(x > 0 for x in row)
                push!(rowinds, rowind)
            end
        end
        rowind += 1
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

#############################################################################################
#############################################################################################
#############################################################################################
############################CORRECTED FUNCTIONS FOR CRITICAL RADIUS##########################

"""
Takes a number and an index position i, and zeros out all digits after the specified position.
The decimal point is not counted towards the position.
"""
function zero_out_after_position(num::Float64, i::Int)
    # Convert number to string
    num_str = string(num)
    
    # Find the decimal point position if it exists
    decimal_pos = findfirst('.', num_str)
    total_n_digits = length(num_str) - 1
    n_digits_whole = length(num_str[1:(decimal_pos-1)])
    n_digits_decimal = length(num_str[(decimal_pos+1):end])
    
    if i == 0
        #the index is 0, don't zero out any elements
        return num
    elseif i == decimal_pos
        elements_after_dot = i - n_digits_whole
        #i is at the decimal point or further down
        return parse(Float64, num_str[1:n_digits_whole]*"."*num_str[(decimal_pos+1):(decimal_pos+elements_after_dot)])
    elseif i < decimal_pos 
        #the last non-zero comes before the decimal position, its easy:
        return parse(Float64, num_str[1:i]*"0"^(n_digits_whole - i))
    else 
        #i is beyond the decimal place
        elements_after_dot = i - n_digits_whole
        return parse(Float64, num_str[1:n_digits_whole]*"."*num_str[(decimal_pos+1):(decimal_pos+elements_after_dot)])
    end
end

"""
Given two numbers, returns the digit position in which the two numbers first differ,
not counting the decimal point.
"""
function first_different_digit_position(num1::Float64, num2::Float64)
    # Convert numbers to strings to compare digit by digit
    str1 = string(num1)
    str2 = string(num2)

    # Remove the decimal point from both strings
    str1_no_dec = replace(str1, "." => "")
    str2_no_dec = replace(str2, "." => "")

    # Get the length of the shorter string without the decimal
    min_length = min(length(str1_no_dec), length(str2_no_dec))

    # Compare digits
    for i in 1:min_length
        if str1_no_dec[i] != str2_no_dec[i]
            return i  # Return the 1-based index of the first differing digit
        end
    end

    # If all compared digits are the same, check if one number has more digits
    if length(str1_no_dec) != length(str2_no_dec)
        return min_length + 1  # Return position after the last matched digit
    end

    return 0  # Return 0 if numbers are identical
end

"""
distinguish_decimal(number::Float64, precision::Float64) -> Float64

Returns the number without the part of the decimal that is below the specified precision.
For example, if the decimal part is 0.23 and the precision is 0.1, it returns
the integer part plus 0.2.
"""
function distinguish_decimal(number::Float64, prec::Float64)
    sign = number < 0 ? -1.0 : 1.0  # Record the sign
    abs_number = abs(number)        # Work with the absolute value

    whole_part = floor(abs_number)   # Get the whole number part
    decimal_part = abs_number - whole_part  # Get the decimal part

    digits_to_keep = first_different_digit_position(decimal_part, decimal_part - prec)
    #round up to number of identical postitions plus one
    truncated = zero_out_after_position(abs_number, digits_to_keep)
    return sign * (truncated)  # Return with the original sign
end

"""
Given a "nperturbations" by "n" matrix, determine if all the rows are positive 
(that is, all have only elements greater than 0), and different from each other 
(that is, each row is unique within a certain tolerance). 
If any row is negative, returns -1
If any row is not unique, returns 0
If neither occurs, returns 1
"""
function check_rows_real_positive_and_unique(matrix::Matrix{ComplexF64}, tolerance::Float64)
    #go through all elements to verify if they are real and positive
    for row in eachrow(matrix)  # Loop through each row
        for element in row  # Loop through each element in the row
            # Check if imaginary part is within tolerance and real part is positive
            if abs(imag(element)) > tolerance || real(element) <= 0 || isnan(real(element))
                return -1  # Return false if any element doesn't meet the criteria
            end
        end
    end
    #if we get here it means all rows are postive and real.
    #now we check that they are unique
    #transform the matrix to a real matrix
    matrix = real.(matrix)
    # Set to track unique rows (rounded for tolerance)
    unique_rows = Set{Vector{Float64}}()

    for row in eachrow(matrix)
        # Create a rounded version of the row for comparison
        row_up_to_precision =  [distinguish_decimal(element, tolerance) for element in row]
        # Check if the rounded row is already in the set
        if row_up_to_precision in unique_rows
            return 0  # Not unique
        end

        # Add the rounded row to the set
        push!(unique_rows, row_up_to_precision)
    end
    #if we get here, we have passed all tests, so we return 1
    return 1
end

"""
gets the first set of perturbations. all must be positive, but different from each other (within tolerance)
"""
function get_first_target(parameters::Tuple, 
                          rho::Float64, 
                          perturbations::Matrix{Float64}, 
                          nperturbations::Int64, 
                          x::Vector{Variable},
                          n::Int64, 
                          tol::Float64)
    perturbed_equilibria = perturbondisc(perturbations, rho, parameters, n, x, false, false, "follow")
    #confirm that indeed they are different and positive
    first_target = select_feasible_equilibria(perturbed_equilibria, 
                                              ones(nperturbations, n))
    is_target_valid = check_rows_real_positive_and_unique(first_target, tol)
    while is_target_valid  != 1
        if is_target_valid == -1
            #some rows are negative or complex, reduce perturbation
            rho = rho/2
        else
            #at least one row is not unique, increase perturbation
            rho = 2*rho
        end
        perturbed_equilibria = perturbondisc(perturbations, rho, parameters, n, x, false, false, "follow")
        is_target_valid = check_rows_real_positive_and_unique(first_target, tol)
        first_target = select_feasible_equilibria(perturbed_equilibria, ones(nperturbations, n))
    end 
    return real.(first_target)
end

"""
once solution is found, run backward checks to confirm it is indeed one
"""
function confirm_solution(rhob::Float64, 
                          perturbations::Matrix{Float64},
                          pars::Tuple,
                          n::Int64,
                          x::Vector{Variable},
                          nperturbations::Int64,
                          target::Matrix{Float64}, 
                          mode::String, 
                          tol::Float64)
    #initialize xcheck   
    xcheck = 1.0e6
    for rho in range(rhob, tol, 10)
        perturbed_equilibria = perturbondisc(perturbations, rho, pars, n, x, true, true, mode)
        #some could be non-feasible, if mode == "follow
        xcheck = get_minimum(perturbed_equilibria, nperturbations, target)
        #println("Backward checking: ", rho, "x = ", xcheck)
        #recalculate if negative to make sure it's not a mistake of the package (sometimes it happens)
        if xcheck < -tol
            equilibria = perturbondisc(perturbations, rho, pars, n, x, true, true, mode)
            xcheck = get_minimum(perturbed_equilibria, nperturbations, target)
        end
        #if at some point x becomes negative or complex again, then another 0 exists
        if xcheck < -tol
            return xcheck
        end
    end
    return xcheck
end 

"""
Perform a disc perturbation of the growth rate in `pars` of radius `rho`, 
and get the equilibrium responses after perturbation. When `interrupt == true`, 
the program will stop whenever there is a perturbation generating any non-feasible solutions.
Equilibria will be stored as Matrix{Complex{Float64}} when only_real_sols is false.
if mode == follows, the system is solved using a parameter homotopy that allows tracking
trajectories of the focal equilibrium
"""
function perturbondisc(perturbations::Matrix{Float64},
                       rho::Float64,
                       parameters::Tuple,
                       n::Int64,
                       x::Vector{Variable},
                       interrupt::Bool,
                       only_real_sols::Bool=true,
                       mode::String="all")
    # Initialize equilibria based on the only_real_sols flag
    equilibria = only_real_sols ? Vector{Matrix{Float64}}() : Vector{Matrix{Complex{Float64}}}()
    nperts = size(perturbations, 1)
    # Project perturbations to current radius
    perturbations_rho = rho * perturbations
    for i in 1:nperts
        pert_i = perturbations_rho[i, :]

        # Solve system (get solutions based on the only_real_sols flag)
        solmat = makeandsolve(x, parameters, n, pert_i, only_real_sols, mode)

        # Count solutions (real, feasible)
        nsols = length(solmat)
        nfeas_sols = length(getfeasrowinds(solmat))

        # Handle the interrupt condition
        if only_real_sols && (nsols == 0 || nfeas_sols == 0) && interrupt
            # End perturbation if desired, whenever feasibility breaks
            return equilibria
        end
        
        # Add equilibria to total
        push!(equilibria, solmat)  # Ensure solmat is stored as Float64
    end
    return equilibria
end

"""
Auxiliary function to find the closest row in a matrix to a given vector.
Both the vector and the matrix can be either Floats or Complex numbers.
"""
function closest_row(vector::Vector{T}, matrix::Matrix{U}) where {T, U}
    # Convert the vector to complex if it's not already
    if T <: Real
        vector = Complex{T}.(vector)
    end

    # Convert the matrix to complex if it's not already
    if U <: Real
        matrix = Complex{U}.(matrix)
    end

    min_distance = Inf
    closest_row_index = 0    
    # Iterate through each row of the matrix
    for row_index in 1:size(matrix, 1)
        current_row = matrix[row_index, :]
        # Calculate the least squares distance (magnitude of the difference)
        distance = sum(abs2.(current_row .- vector))

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
function positive_row(target_vector::AbstractVector, matrix::AbstractMatrix)

    positive_rows_indices = getfeasrowinds(matrix)

    # Select the appropriate row based on the positive rows found
    if length(positive_rows_indices) == 1
        return positive_rows_indices[1]  # Return the only positive row found
    elseif length(positive_rows_indices) > 1
        # If multiple positive rows are found, form a matrix of those rows
        positive_matrix = matrix[positive_rows_indices, :]

        # Use the closest_row function to find the closest row to the target vector
        closest_row_index = closest_row(target_vector, positive_matrix)

        return positive_rows_indices[closest_row_index]  # Return the index from the original matrix        
    elseif size(matrix, 1) == 0
        #matrix is empty, return nothing
        return nothing
    else
        # If no positive rows found, use the original matrix to find the closest row
        return closest_row(target_vector, matrix)  # Return closest row index from original matrix
    end
end

"""
Given a collection of equilibra and a target, select equilibria as follows:
Aim to return all positive rows.
If no positive rows exist for a given perturbation, return the closest row to th
target equilibria. 
"""
function select_feasible_equilibria(array_of_matrices::Vector{<:AbstractMatrix},
                                    target_matrix::AbstractMatrix)
    k = length(array_of_matrices)   # Number of matrices

    # Determine the type of the first matrix to initialize matched_matrix accordingly
    first_matrix = array_of_matrices[1]
    matched_matrix = similar(first_matrix, k, size(target_matrix, 2))  # Initialize matched matrix

    for i in 1:k
        # Get the ith row of the target matrix
        target_row = target_matrix[i, :]
        matrix_i = array_of_matrices[i]

        if size(matrix_i, 1) == 0
            # Handle empty matrix case
            matched_matrix[i, :] .= NaN  # or zeros
            continue  # Skip further processing for this iteration
        end

        # Find the closest positive row
        row_index = positive_row(target_row, matrix_i)
        row_i = matrix_i[row_index, :]
        if row_index == 0 || (row_i isa Vector{Complex{Float64}})
            row_index = closest_row(target_row, matrix_i)  # Fallback to follow mode
        end
        matched_matrix[i, :] = matrix_i[row_index, :]
    end
    return matched_matrix
end

"""
given the array of matrices with equilibria from all perturbations, select the minimum of them
after some filtering based on the mode perturbations and result of perturbations
This function is only called inside findmaxperturbation, which assumes interrupt == true always
therefore, the selected_equilibria should always contain positive rows if the mode of selection is
positive. Otherwise the perturbation protocol should have not finished, and perturbed_equilibria should
be shorter than the number of perturbations.
"""
function get_minimum(perturbed_equilibria::Vector{Matrix{Float64}},
                     nperturbations::Int64,
                     target_equilibria::Matrix{Float64}
                     )
    if length(perturbed_equilibria) < nperturbations
        #the perturbation protocol stopped early (non-feasibility encountered)
        xmin = -1.0
    else 
        #all perturbations yielded feasible equilibria
        selected_equilibria = select_feasible_equilibria(perturbed_equilibria,
                                                         target_equilibria)
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
                             n::Int64,
                             x::Vector{Variable},
                             target_equilibria::Matrix{Float64},
                             mode::String,
                             tol::Float64=1e-9
                             )
    #find the middle point of current interval
    global rhob = bisect(rho1, rho2)
    #perform disc perturbation of radius rb
    while abs(rho1 - rho2) > tol
        perturbed_equilibria = perturbondisc(perturbations, rhob, parameters, n, x, true, true, mode)
        #if perturbation protocol stopped early, there were non-feasible solutions
        if length(perturbed_equilibria) < nperturbations
            #the perturbation protocol stopped early (non-feasibility encountered)
            xmin = -1.0
        else 
            #otherwise there are at least one feasible equilibria per perturbation
            selected_equilibria = select_feasible_equilibria(perturbed_equilibria, target_equilibria)
            #this minimum could be negative after selection if mode=="follow"
            xmin = minimum(selected_equilibria)
        end
        #println("[", rho1, ",", rho2, "] with xmin = ", xmin, " at rhob = ", rhob)
        #redefine limits
        rho1, rho2 = getlims(rho1, rho2, rhob, xmin)
        #if solution is found, check that no other solutions exist for smaller rs
        if abs(rho1 - rho2) < tol
            xmin = confirm_solution(rhob, perturbations, parameters, n, x, 
                                    nperturbations, target_equilibria, mode, tol)
            if xmin < -tol
                rho1 = 0.0
                rho2 = rhob
            end
        end
        return findmaxperturbation(rho1, rho2, perturbations, nperturbations, parameters,
                                   n, x, target_equilibria, mode, tol)
    end
    return rhob
end

"""
function to count the average number of positve rows in each matrix of equilibria arising from 
a perturbation.
"""
function average_number_positive_rows(array_of_matrices::Vector{<:AbstractMatrix})
    k = length(array_of_matrices)   # Number of matrices
    # Initialize vector of number of positive rows in each matrix
    n_positive_rows = []  
    for i in 1:k #loop through each matrix
        matrix_i = array_of_matrices[i]
        # Iterate through each row of the matrix
        feasible_inds = getfeasrowinds(matrix_i)
        push!(n_positive_rows, length(feasible_inds))
    end
    return mean(n_positive_rows)
end

function testrmax(perturbations, parameters, n, x, interrupt, rmax, mode)
    # Define the perturbation vector
    println("Rmax is: ", rmax)
    perturbation_vector = collect((.85 * rmax):(0.05 * rmax):(1.15 * rmax))

    # Initialize an array to hold the resulting matrices
    results = []

    for r in perturbation_vector
        # Call perturbondisc for each perturbation
        result = perturbondisc(perturbations, r, parameters, n, x, false, mode)
        flatresult = vcat(result...)
        push!(results, vcat(result...))
    end
    # Flatten the results into a single matrix
    flattened_results = vcat(results...)

    # Plot the points
    display(plot(flattened_results[:, 1], 
         flattened_results[:, 2], 
         seriestype = :scatter, 
         label = "Perturbation Points", 
         legend = false,
         xlimits = (0, maximum(flattened_results[:,1])),
         ylimits = (0, maximum(flattened_results[:,2]))))
end

#I should write a better test of rmax. one that increases r right above the tolerance and checks that indeed 
#complex or negative solutions appear.