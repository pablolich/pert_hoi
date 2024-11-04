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

function test_getfeasrowinds()
    # Test with a real matrix
    real_matrix = [
        1.0 2.0;  # Feasible
        -1.0 3.0; # Not feasible
        0.5 0.1;  # Feasible
        0.0 2.0   # Not feasible
    ]
    expected_real_indices = [1, 3]  # Only rows 1 and 3 are feasible
    actual_real_indices = getfeasrowinds(real_matrix)

    println("Testing with real matrix...")
    println("Expected: ", expected_real_indices)
    println("Actual: ", actual_real_indices)
    @assert actual_real_indices == expected_real_indices

    # Test with a complex matrix
    complex_matrix = [
        1 + 1im 2 + 2im;      # Not feasible (imaginary parts > tol)
        0.5 + 0im 0.1 + 0im;  # Feasible (real parts > 0)
        1 - 1e-11im 2 + 0im;  # Feasible (imaginary parts < tol)
        -1 + 0im 3 + 1im      # Not feasible (negative real part)
    ]
    expected_complex_indices = [2, 3]  # Rows 2 and 3 are feasible
    actual_complex_indices = getfeasrowinds(complex_matrix)

    println("Testing with complex matrix...")
    println("Expected: ", expected_complex_indices)
    println("Actual: ", actual_complex_indices)
    @assert actual_complex_indices == expected_complex_indices

    println("All tests passed!")
end

function test_makeandsolve()
    n = 2  # Define n first
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n

    result = makeandsolve(x, parameters, n)
    # Check if the result is a Vector of Matrix{Float64}
    if isa(result, Matrix{Float64})
        println("Test makeandsolve: Passed - Output is a Vector of Matrix{Float64}")
    else
        println("Test makeandsolve: Failed - Output type is ", typeof(result))
    end
    #now test the output when real = false
    matrix_solutions = makeandsolve(x, parameters, n, false)
    println("matrix sols with complex ", matrix_solutions)
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

    result = perturbondisc(perturbations, rho, parameters, n, x, interrupt)

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

function test_closest_row()
    n = 2  # Define n first
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    #now test the output when real = false
    matrix_solutions = makeandsolve(x, parameters, n, false)
    vector = [1.0, 1.0]

    result = closest_row(vector, matrix_solutions)

    println("Test closest_row: Inputed vector is ", vector)
    println("Inputed matrix is ", matrix_solutions)
    println("Closest row is ", matrix_solutions[result,:])
end

function test_positive_row()
    n = 2
    # Test cases
    test_cases = [
        # Test case 1: No positive rows
        ([1.0, 2.0], [0.0 0.0; 
                      -1.0 -2.0; 
                      1.1 2.1
                      -3.0 -4.0;
                      1.2 2.2], 3),  # Closest row should be the thirdx
        # Test case 2: One positive row
        ([1.0, 2.0], [0.0 0.0; 1.0 1.0], 2),  # Should return the index of the positive row
        # Test case 3: Multiple positive rows
        ([1.0, 2.0], [1.0 2.0; 3.0 4.0; 2.0 2.0], 1),  # Closest to [1,2] should be the first row
        # Test case 4: All negative elements
        ([1.0, 1.0], [-1.0 -2.0; -3.0 -4.0], 1),  # Closest should be the first (since all are negative)
        # Test case 5: Empty matrix
        ([1.0, 1.0],Matrix{Float64}(undef, 0, n), nothing)  # Should handle gracefully
    ]

    for (target_vector, matrix, expected) in test_cases
        result = positive_row(target_vector, matrix)
        if expected == nothing
            println("Test with matrix $matrix: Expected nothing, got $result - ", result == nothing ? "Passed" : "Failed")
        else
            println("Test with matrix $matrix: Expected index $expected, got $result - ", result == expected ? "Passed" : "Failed")
        end
    end
end

function test_select_equilibria()
    n = 2  # Define n first
    array_of_matrices = [
        [1.0 2.0; 0.1 4.0; 1 0.1], #perturbation 1
        Matrix{Float64}(undef, 0, n), #perturbation 2
        [0.1 0.1; 7.0 8.0; 1.1 .1], #perturbation 3
    ]
    target_matrix = [1.0 2.0; 5.0 6.0; 1 -1]
    mode = "follow"

    result = select_equilibria(array_of_matrices, target_matrix, mode)
    println("Test select_equilibria (mode=follow): Expected selected equilibria, got: ", result, " - Passed")

    mode = "positive"
    result = select_equilibria(array_of_matrices, target_matrix, mode)
    print("Test select_equilibria (mode=positive): Expected matrix, got ", result, "- Passed")

    # New test for complex matrices with non-zero imaginary parts
    complex_array_of_matrices = [
        [Complex{Float64}(1.0, 1.0) Complex{Float64}(2.0, -1.0);
         Complex{Float64}(0.1, 0.5) Complex{Float64}(4.0, 0.0);
         Complex{Float64}(1.0, -0.2) Complex{Float64}(0.1, 3.0)], # Complex perturbation 1
        Matrix{Complex{Float64}}(undef, 0, n), # Complex perturbation 2
        [Complex{Float64}(0.1, 0.1) Complex{Float64}(0.1, -0.3);
         Complex{Float64}(7.0, 0.4) Complex{Float64}(8.0, 1.0);
         Complex{Float64}(1.1, -0.5) Complex{Float64}(0.1, 0.2)], # Complex perturbation 3
    ]
    complex_target_matrix = [Complex{Float64}(0.1, 0.4) Complex{Float64}(4.3, -0.1);
                             Complex{Float64}(5.0, 0.5) Complex{Float64}(6.0, 1.0);
                             Complex{Float64}(1.0, -1.0) Complex{Float64}(0.0, 0.0)]
    complex_mode = "follow"

    complex_result = select_equilibria(complex_array_of_matrices, complex_target_matrix, complex_mode)
    println("Test select_equilibria (complex, mode=follow): Expected selected equilibria, got: ", complex_result)
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
    tol = 1e-9

    try
        result = get_first_target(parameters, rho, perturbations, nperturbations, x, n, tol)
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
    rho2 = 10.0
    perturbations = points_hypersphere(2, 1.0, 10)
    nperturbations = size(perturbations, 1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    target = get_first_target(parameters, 1e-8, perturbations, nperturbations, x, n, tol)
    mode = "follow"
    # Test with initial rho1
    rhomax, equilibria = findmaxperturbation(rho1, rho2, perturbations, nperturbations, parameters, n, x, target, mode, tol)
    #display(plotperturbations(equilibria))
    if isa(rhomax, Float64) && isa(equilibria, Matrix{Float64}) && size(equilibria, 1) == nperturbations
        println("Test findmaxperturbation: Passed - Expected float and matrix, got: ", typeof(rhomax), " ", typeof(equilibria), size(equilibria))
    else
        println("Test findmaxperturbation: Failed - Expected a float and matrix, got: ", typeof(rhomax), typeof(equilibria), size(equilibria))
    end

    # Test for slightly increased rhomax
    rho_test = rhomax + 0.0001

    # Test perturbondisc with increased rho
    interrupt = true
    result_interrupt = perturbondisc(perturbations, rho_test, parameters, n, x, interrupt)
    println(result_interrupt)
    selected_equilibria = select_equilibria(result_interrupt, equilibria, "follow")
    println(selected_equilibria)
    #display(plotperturbations(selected_equilibria))
    println(selected_equilibria)
    if size(selected_equilibria, 1) < nperturbations
        println("Test perturbondisc (interrupt=true): Passed - Output has fewer rows than nperturbations: ", size(result_interrupt, 1))
    else
        println("Test perturbondisc (interrupt=true): Failed - Expected fewer rows than nperturbations, got: ", size(result_interrupt, 1))
    end

#    interrupt = false
#    result_no_interrupt = perturbondisc(perturbations, rho_test, parameters, x, interrupt)

#    # Check if there are any negative rows
#    negative_rows = count(row -> any(row .< 0), eachrow(result_no_interrupt))
#    if negative_rows > 0
#        println("Test perturbondisc (interrupt=false): Passed - Found negative rows: ", negative_rows)
#    else
#        println("Test perturbondisc (interrupt=false): Failed - No negative rows found.")
#    end
end

function test_testrmax()
    n = 2  # Define n first
    rng = MersenneTwister(1)

    tol = 1e-9
    rho1 = tol
    rho2 = 10.0
    perturbations = points_hypersphere(2, 1.0, 100)
    nperturbations = size(perturbations, 1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    mode="follow"

    target = get_first_target(parameters, 1e-3, perturbations, nperturbations, x, n, tol)
    #find rmax 
    rmax, equilibria = findmaxperturbation(tol, 10.0, perturbations, nperturbations, parameters, n, x, target, mode, tol)
    testrmax(perturbations, parameters, n, x, false, rmax)

end

function test_average_number_positive_rows()
    n = 2  # Define n first
    perturbations = points_hypersphere(2, 1.0, 100)
    rho = 0.004
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    interrupt = true

    response_perts = perturbondisc(perturbations, rho, parameters, n, x, interrupt)

    result = average_number_positive_rows(response_perts)
    println("Test average_number_positive_rows: Expected average number of positive rows, got: ", result, " - Passed")

    #now test with complex entries
    response_perts = perturbondisc(perturbations, rho, parameters, n, x, false, false)
    result = average_number_positive_rows(response_perts)
    println("Test average_number_positive_rows with complex entries: Expected average number of positive rows, got: ", result, " - Passed")

end

function test_proportion_of_positive_rows()
    n = 2  # Define n first
    perturbations = points_hypersphere(2, 1.0, 100)
    rho = 0.004
    rng = MersenneTwister(1)
    parameters_nested = (0.5, sampleparameters(n, rng, 1))  # Adjusted values in the tuple
    # Flatten tuple
    parameters = Tuple(Iterators.flatten(parameters_nested))
    @var x[1:n]  # Define x as a function of n
    matrix = perturbondisc(perturbations, rho, parameters, n, x, false, false)

    result = proportion_of_positive_rows(real(matrix))
    println("Test proportion_of_positive_rows: Expected proportion of positive rows, got: ", result, " - Passed")

    result_alternative = getfeasrowinds(matrix)

    @assert result == length(result_alternative)/size(matrix, 1)
end

# Call test functions
# println("Test fzero_out_after_position()")
# test_zero_out_after_position()

# println()
# println("Test first_different_digit_position()")
# test_first_different_digit_position()

# println()
# println("Test distinguish_decimal()")
# test_distinguish_decimal()

# println()
# println("Test check_rows_positive_and_unique()")
# test_check_rows_positive_and_unique()

# println()
# println("Test getfeasrowinds()")
# test_getfeasrowinds()

# println()
# println("Test makeandsolve()")
# test_makeandsolve()

# println()
# println("Test sampleparameters()")
# test_sampleparameters()

# println()
# println("Test perturbondisc()")
# test_perturbondisc()

# println()
# println("Test closest_row()")
# test_closest_row()

# println()
# println("Test postive_row)")
# test_positive_row()

println()
println("Test select_equilibria()")
test_select_equilibria()

# println()
# println("Test get_minimum()")
# test_get_minimum()

# println()
# println("Test confirm_solution()")
# test_confirm_solution()

# println()
# println("Test get_first_target")
# test_get_first_target()


# println()
# println("Test findmaxperturbation()")
# test_findmaxperturbation()

# println()
# println("Test testrmax()")
# test_testrmax() 

# println()
# println("Test average_number_positive_rows()")
# test_average_number_positive_rows()

# println()
# println("Test proportion_of_positive_rows()")
# test_proportion_of_positive_rows()