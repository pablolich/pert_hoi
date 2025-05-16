using DelimitedFiles, LinearAlgebra, Random
using GLMakie

# Generate n random points on unit sphere in 3D
function random_unit_sphere_points(n::Int)
    points = randn(n, 3)
    points ./= sqrt.(sum(points.^2, dims=2))  # normalize each row
    return points
end

# Load and reproject optimized points
function load_and_reproject_points(filename::String)
    pts = readdlm(filename)
    pts ./= sqrt.(sum(pts.^2, dims=2))  # normalize each point
    return pts
end

# Generate a sphere surface as x, y, z matrices
function generate_sphere_surface(resolution=50)
    θ = range(0, π, length=resolution)
    ϕ = range(0, 2π, length=resolution)
    x = [sin(t)*cos(p) for t in θ, p in ϕ]
    y = [sin(t)*sin(p) for t in θ, p in ϕ]
    z = [cos(t) for t in θ, _ in ϕ]
    return x, y, z
end

# Plot random and optimized points on unit sphere
function plot_points_comparison(random_pts, optimized_pts)
    fig = Figure(size=(1000, 500))

    x_sphere, y_sphere, z_sphere = generate_sphere_surface()

    ax1 = LScene(fig[1, 1], show_axis=false)
    surface!(ax1, x_sphere, y_sphere, z_sphere, color=:gray, shading=false, transparency=true, alpha=0.3)
    scatter!(ax1, random_pts[:,1], random_pts[:,2], random_pts[:,3], markersize=5, color=:blue)
    Label(fig[0, 1], "Random Points", tellwidth=false)

    ax2 = LScene(fig[1, 2], show_axis=false)
    surface!(ax2, x_sphere, y_sphere, z_sphere, color=:gray, shading=false, transparency=true, alpha=0.3)
    scatter!(ax2, optimized_pts[:,1], optimized_pts[:,2], optimized_pts[:,3], markersize=5, color=:blue)
    Label(fig[0, 2], "Optimized Points", tellwidth=false)

    return fig
end

# Run
random_pts = random_unit_sphere_points(500)
optimized_pts = load_and_reproject_points("../../data/thompson_perturbations/optimized_n500_d3_cost_115802.750913.txt")
fig = plot_points_comparison(random_pts, optimized_pts)
save("sphere_points_comparison.png", fig)
