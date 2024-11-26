#sample parameters n = 2
#pick one alpha
#get the real feasibility domain
#get critical radius with many points (1e4)
#get critical radius with a fraction of points (1e2)
#plot the r's coming from the real feasibility domain
#plot a circle with the radius equal the one obtained for 1e4 (i.e., inf)
#plot a the points used to find critical radius of 1e2 points
#if desired, plot also the x corresponding to this cases

"""
find middle point between two points
"""
function bisect(x0::Float64, x1::Float64)
    return 1/2*(x0 + x1)
end

"""
get appropriate limits based on sign of the middle point
"""
function getlims(x0::Float64, x1::Float64, xb::Float64, xmin::Float64)
    if xmin > 0
        x0 = xb
    else
        x1 = xb
    end
    return x0, x1
end

"""
get 
"""
function get_angles(num_points::Int64)
    points = zeros(num_points, 2)  # Initialize matrix to store points
    for i in 1:num_points
        theta = 2π * (i - 1) / num_points  # Calculate angle for current point
        points[i, 1] = cos(theta)  # Calculate x coordinate
        points[i, 2] = sin(theta)  # Calculate y coordinate
    end
    return points
end

"""
get the actual critical perturbation boundary for which feasibility is lost
"""
function get_critical_domain(parameters::Tuple, )
    #pick a direction and increase r until feasibility is lost
    #if feasibility is not lost, then pick another (random) direction and do the same thing
    #to see determine when is feasibility lost, use bisection method.
    #in particular, find the r that makes at least one component of x 0 or complex
    #once that point is found, rotate the r slightly and perform the search again.
    #the limits of the search should be a 15% increase on both directions from the previous search
    return 
end
