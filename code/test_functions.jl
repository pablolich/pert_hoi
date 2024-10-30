include("functions.jl")

# Test functions
function test_zero_out_after_position()
    test_cases = [
        (123.4567, 5, 123.4500),   # Zero out after 5th position
        (5.678, 1, 5.000),          # Zero out after 1st position
        (9.87654, 3, 9.87000),      # Zero out after 3rd position
        (2.5, 2, 2.500),            # No change as i is at the end
        (100.999, 4, 100.9000),      # Zero out after 4th position
        (100.333, 1, 100.0)
    ]

    for (num, i, expected) in test_cases
        result = zero_out_after_position(num, i)
        println("Zeroing out after position $i for $num => Result: $result (Expected: $expected) - ", 
                result == expected ? "Passed" : "Failed")
    end
end

function test_first_different_digit_position()
    test_cases = [
        (5.23, 5.24, 3),      # First difference at position 3
        (5.567, 5.56, 4),     # First difference at position 4
        (10.0, 10.1, 3),      # First difference at position 3
        (3.15, 3.15, nothing), # Identical numbers
        (2.0001, 2.001, 4)     # First difference at position 4
    ]

    for (num1, num2, expected) in test_cases
        result = first_different_digit_position(num1, num2)
        println("Comparing $num1 and $num2 => Result: $result (Expected: $expected) - ", 
                result == expected ? "Passed" : "Failed")
    end
end

function test_distinguish_decimal()
    test_cases = [
        (5.23, 0.1, 5.2),
        (5.567, 0.05, 5.56),
        (5.02, 0.05, 5.0),
        (5.09, 0.1, 5.0),
        (3.15, 0.1, 3.1),
        (3.567, 0.1, 3.5),
        (-2.03, 0.05, -2.0),
        (-5.23, 0.1, -5.2),
        (-5.567, 0.05, -5.56),
        (10.0, 0.1, 10.0),   # Testing whole number
        (10.005, 0.01, 10.0)  # Decimal part below precision
    ]

    for (number, precision, expected) in test_cases
        result = distinguish_decimal(number, precision)
        println("Input: $number with precision $precision => Result: $result (Expected: $expected) - ", 
                result == expected ? "Passed" : "Failed")
    end
end

function test_check_rows_positive_and_unique()
    matrix1 = [1.0 2.0; 3.0 4.0; 5.0 6.0]
    matrix2 = [1.0 2.0; 1.011 2.019; 5.0 6.0]
    matrix3 = [1.0 2.0; -1.0 4.0; 5.0 6.0]

    println("Test check_rows_positive_and_unique (should be true): ", 
            check_rows_positive_and_unique(matrix1, 0.01) ? "Passed" : "Failed")
    println("Test check_rows_positive_and_unique (should be true): ", 
            check_rows_positive_and_unique(matrix2, 0.01) ? "Passed" : "Failed")
    println("Test check_rows_positive_and_unique (should be false): ", 
            !check_rows_positive_and_unique(matrix3, 0.01) ? "Passed" : "Failed")
end

function test_makeandsolve()
    n = 2  # Define n first
    parameters = (2.0, [1.5, 1.5], [0.1 0.2; 0.3 0.4], rand(Float64, n, n, n))  # Adjusted values in the tuple
    rho = 0.1
    @var x[1:n]  # Define x as a function of n

    result = makeandsolve(x, parameters, n)
    # Check if the result is a Vector of Matrix{Float64}
    if isa(result, Matrix{Float64})
        println("Test makeandsolve: Passed - Output is a Vector of Matrix{Float64}")
    else
        println("Test makeandsolve: Failed - Output type is ", typeof(result))
    end
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

function test_sampleparameters()
    n = 2  # Define n first
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n

    matrix_solutions = makeandsolve(x, parameters, n)
    result = has_row_all_ones(matrix_solutions)
    if result
        println("Test sampleparameters: Passed -- Matrix of solutions contains at least one row with all ones")
    else
        println("Test sampleparameters: Failed -- Matrix of solutions does not contain at least one row with all ones")
    end
end

function test_perturbondisc()
    n = 2  # Define n first
    perturbations = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    rho = 0.1
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    interrupt = true

    result = perturbondisc(perturbations, rho, parameters, x, interrupt)

    # Check if the result is a Vector of Matrix{Float64}
    if isa(result, Vector{Matrix{Float64}})
        println("Test perturbondisc: Passed - Output is a Vector of Matrix{Float64}")
        
        # Additional check for interrupt = true
        if interrupt
            all_positive_rows = all(matrix -> any(all(row .> 0) for row in eachrow(matrix)), result)
            if all_positive_rows
                println("Test perturbondisc: Passed - Every matrix contains at least one row with all positive elements.")
            else
                println("Test perturbondisc: Failed - Not every matrix contains a row with all positive elements.")
            end
        end
    else
        println("Test perturbondisc: Failed - Output type is ", typeof(result))
    end
end

function test_select_equilibria()
    n = 2  # Define n first
    array_of_matrices = [
        [1.0 2.0; -1.0 4.0; 3.0 4.0],
        Matrix{Float64}(undef, 0, n),
        [5.0 6.0; 7.0 8.0; 9.0 10.0]
    ]
    target_matrix = [1.0 2.0; 5.0 6.0; 1 -1]
    mode = "follow"

    result = select_equilibria(array_of_matrices, target_matrix, mode)
    println("Test select_equilibria: Expected selected equilibria, got: ", result, " - Passed")
end

function test_get_minimum()
    n = 2  # Define n first
    all_equilibria = [
        [1.0 2.0; 3.0 4.0; 5.0 6.0],
        [7.0 8.0; 9.0 10.0; 11.0 12.0]
    ]
    nperturbations = 2
    target_equilibria = [1.0 2.0; 7.0 8.0]
    mode = "follow"

    result = get_minimum(all_equilibria, nperturbations, target_equilibria, mode)
    println("Test get_minimum: Expected minimum value, got: ", result, " - Passed")
end

function test_confirm_solution()
    n = 2  # Define n first
    perturbations = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    rhob = 0.1
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    nperturbations = 3
    target = [1.0 2.0; 3 4; 5 7]
    mode = "follow"
    tol = 1e-5
    try
        result = confirm_solution(rhob, perturbations, parameters, n, x, nperturbations, target, mode, tol)
        println("Test confirm_solution: Expected a float64, got: ", result, " - Passed")
    catch e
        println("Test confirm_solution: Expected no error, caught exception: ", e, " - Failed")
    end
end

function test_get_first_target()
    # Mock data for test
    n = 2  # Define n first
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))    
    rho = 0.1
    perturbations = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    nperturbations = size(perturbations, 1)
    @var x[1:n]  # Define x as a function of n

    try
        result = get_first_target(parameters, rho, perturbations, nperturbations, x, n)
        println("Test get_first_target: Expected a matrix of equilibria, got: ", result, " - Passed")
    catch e
        println("Test get_first_target: Expected no error, caught exception: ", e, " - Failed")
    end
end

function test_findmaxperturbation()
   
    n = 2  # Define n first
    rng = MersenneTwister(1)

    tol = 1e-9
    rho1 = tol
    rho2 = 10
    perturbations = [1.0 1.0; 2.0 2.0; 3.0 3.0]
    nperturbations = 3
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    target = [1.0 2.0; 3 4; 5 7]
    mode = "follow"

    try
        result = findmaxperturbation(rho1, rho2, perturbations, nperturbations, parameters, x,target, mode,tol)
        println("Test findmaxperturbation: Expected a float64, got: ", result, " - Passed")
    catch e
        println("Test findmaxperturbation: Expected no error, caught exception: ", e, " - Failed")
    end
end

function test_closest_row()
    n = 2  # Define n first
    vector = [1.0, 2.0]
    matrix = [1.0 2.0; 3.0 4.0; 5.0 6.0]

    result = closest_row(vector, matrix)
    println("Test closest_row: Expected closest row index, got: ", result, " - Passed")
end

function test_positive_row()
    n = 2  # Define n first
    target_vector = [1.0, 2.0]
    matrix = [1.0 2.0; -1.0 4.0; 3.0 4.0]

    result = positive_row(target_vector, matrix)
    println("Test positive_row: Expected positive row index, got: ", result, " - Passed")
end

function test_average_number_positive_rows()
    n = 2  # Define n first
    array_of_matrices = [
        [1.0 2.0; -1.0 4.0; 3.0 4.0],
        [5.0 6.0; 7.0 8.0; -9.0 10.0]
    ]

    result = average_number_positive_rows(array_of_matrices)
    println("Test average_number_positive_rows: Expected average number of positive rows, got: ", result, " - Passed")
end

function test_proportion_of_positive_rows()
    n = 2  # Define n first
    matrix = [1.0 2.0; -1.0 4.0; 3.0 4.0]

    result = proportion_of_positive_rows(matrix)
    println("Test proportion_of_positive_rows: Expected proportion of positive rows, got: ", result, " - Passed")
end

# Call test functions
println("Test fzero_out_after_position()")
test_zero_out_after_position()

println()
println("Test first_different_digit_position()")
test_first_different_digit_position()

println()
println("Test distinguish_decimal()")
test_distinguish_decimal()

println()
println("Test check_rows_positive_and_unique()")
test_check_rows_positive_and_unique()

println()
println("Test makeandsolve()")
test_makeandsolve()

println()
println("Test sampleparameters()")
test_sampleparameters()

println()
println("Test perturbondisc()")
test_perturbondisc()

println()
println("Test select_equilibria()")
test_select_equilibria()

println()
println("Test get_minimum()")
test_get_minimum()

println()
println("Test confirm_solution()")
test_confirm_solution()

println()
println("Test get_first_target")
test_get_first_target()


# println()
# println("Test findmaxperturbation()")
# test_findmaxperturbation()


# test_closest_row()
# test_positive_row()
# test_average_number_positive_rows()
# test_proportion_of_positive_rows()