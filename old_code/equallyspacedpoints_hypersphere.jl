using LinearAlgebra
using Random
using Plots

function spherical_to_cartesian(rho, angles)
    n = length(angles)
    coords = Vector{Float64}(undef, n)
    
    coords[1] = rho * cos(angles[1])
    
    for i in 2:n
        coords[i] = rho * sin(angles[i-1]) * cos(angles[i])
    end
    
    coords[n] = rho * prod(sin.(angles[1:n-1]))
    
    return coords
end

function generate_hypersphere_points(n, rho, N)
    points = Matrix{Float64}(undef, N, n)
    
    phi = pi * (3 - sqrt(5))  # Golden angle in radians
    
    for k in 1:N
        index = k - 0.5
        theta = acos(1 - 2 * index / N)
        phi_angle = phi * (k - 1)
        
        angles = [theta] # First angle is theta
        for i in 1:(n - 1)
            push!(angles, phi_angle * i)
        end
        
        points[k, :] = spherical_to_cartesian(rho, angles)
    end
    
    return points
end

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
calculate euclidean distance between all pairs of points
"""
function all_distances(points)
    npoints = size(points, 1)
    min_distances = []
    for i in 1:npoints
        pointi = points[i,:]
        distancesiall = []
        for j in 1:npoints
            if j == i
                continue
            else
                pointj = points[j,:]
                #calculate distance
                distij = norm(pointi - pointj)
                push!(distancesiall, distij)
            end
        end
        #get minimum distance
        mindisti = minimum(distancesiall)
        push!(min_distances, mindisti)
    end
end


# Example usage
n = 3 # Dimension
rho = 1.0  # Radius
N = 500  # Number of points

points = generate_hypersphere_points(n, rho, N)
points_rand = points_hypersphere(n, rho, N)
min_distances_rand = all_distances(points_rand)
min_distances = all_distances(points)

println(min_distances_rand)
p1 = scatter(points[:,1], points[:,2], points[:, 3], markersize = 2)
p2 = scatter(points_rand[:,1], points_rand[:,2], points_rand[:,3], markersize=2)
p3 = histogram(min_distances)
p4 = histogram(min_distances_rand)
plot(p1, p2, p3, p3)