using LinearAlgebra, Random
using Distributions
using CairoMakie

# Constraint matrix builder
documentation = """
    build_C(ℓ, m, n)

Construct the 2n×D matrix C enforcing:
 1) For each i=1:n, sum over A[i, j] plus r[i] == 0
 2) For each i=1:n, sum over B[:,:,i] plus r[i] == 0
Where D = n + m*n + ℓ*m*n.
"""
function build_C(ℓ::Int, m::Int, n::Int)
    D = n + m*n + ℓ*m*n
    C = zeros(Float64, 2n, D)

    # 1) sum_j A[i,j] + r[i] == 0  (row sums of A)
    for i in 1:n
        C[i, i] = 1.0
        for j in 1:n
            col = n + (j-1)*m + i
            C[i, col] = 1.0
        end
    end

    # 2) sum_{j,k} B[j,k,i] + r[i] == 0
    offsetB = n + m*n
    for i in 1:n
        row = n + i
        C[row, i] = 1.0
        for k in 1:m, j in 1:ℓ
            lin = j + (k-1)*ℓ + (i-1)*ℓ*m
            C[row, offsetB + lin] = 1.0
        end
    end

    return C
end

# Modified hit-and-run sampler using a truncated Normal(0,1) step
"""
    hit_and_run_sample(C, L, U, x0, ℓ, m, n)

Perform one hit-and-run move sampling λ from a truncated Normal(0,1) on the feasible chord.
Returns (r, A, B)
"""
function hit_and_run_sample(
    C::AbstractMatrix{<:Float64},
    L::AbstractVector{<:Float64},
    U::AbstractVector{<:Float64},
    x0::AbstractVector{<:Float64},
    ℓ::Int, m::Int, n::Int;
    max_tries::Int = 100
)
    # Precompute nullspace basis
    F     = svd(C; full=true)
    tol   = maximum(size(C)) * eps(eltype(C)) * maximum(F.S)
    rankC = count(>(tol), F.S)
    Q     = F.V[:, rankC+1:end]
    D     = length(x0)

    for attempt in 1:max_tries
        # random direction in ker(C)
        u = Q * randn(size(Q,2))
        nu = norm(u)
        if nu < 1e-12
            continue
        end
        u ./= nu

        # find feasible interval
        λ_min, λ_max = -Inf, Inf
        for j in 1:D
            if abs(u[j]) > 1e-12
                a = (L[j] - x0[j]) / u[j]
                b = (U[j] - x0[j]) / u[j]
                lo, hi = min(a,b), max(a,b)
                λ_min = max(λ_min, lo)
                λ_max = min(λ_max, hi)
            end
        end
        if λ_max < λ_min || !isfinite(λ_min) || !isfinite(λ_max)
            continue
        end

        # sample λ from truncated Normal(0,1)
        λ = rand() * (λ_max - λ_min) + λ_min
        x = x0 .+ λ .* u

        # split back into (r, A, B)
        r =     x[1                 : n]
        A = reshape(x[n+1             : n + m*n], m, n)
        B = reshape(x[n + m*n + 1     : end],       ℓ, m, n)

        return r, A, B
    end
    error("Failed to find a valid hit-and-run step after \$max_tries tries.")
end

"""
    sample_and_plot_histograms(ℓ, m, n, N; lower=-10.0, upper=10.0, output_file="histograms.png")

Run N hit-and-run samples under the constraints
  C*x == 0, lower ≤ x ≤ upper
for dimensions:
  r ∈ ℝⁿ, A ∈ ℝ^{m×n}, B ∈ ℝ^{ℓ×m×n},
then plot histograms of r, A, and B entries side-by-side.
"""
function sample_and_plot_histograms(
    ℓ::Int, m::Int, n::Int, N::Int;
    lower::Float64 = -10.0,
    upper::Float64 = 10.0,
    output_file::String = "histograms.png"
)
    D = n + m*n + ℓ*m*n
    C = build_C(ℓ, m, n)
    L = fill(lower, D)
    U = fill(upper, D)

    # initial point
    r0 = zeros(n)
    A0 = zeros(m, n)
    B0 = zeros(ℓ, m, n)
    x_current = vcat(r0, vec(A0), vec(B0))

    rs = zeros(n, N)
    As = zeros(m*n, N)
    Bs = zeros(ℓ*m*n, N)

    for i in 1:N
        r, A, B = hit_and_run_sample(C, L, U, x_current, ℓ, m, n)
        rs[:, i] = r
        As[:, i] = vec(A)
        Bs[:, i] = vec(B)
        x_current = vcat(r, vec(A), vec(B))
    end

    fig = Figure(resolution = (1200, 400))
    ax1 = Axis(fig[1, 1], title = "r entries", xlabel = "Value", ylabel = "Count")
    hist!(ax1, vec(rs), bins = 50)
    ax2 = Axis(fig[1, 2], title = "A entries", xlabel = "Value", ylabel = "Count")
    hist!(ax2, vec(As), bins = 50)
    ax3 = Axis(fig[1, 3], title = "B entries", xlabel = "Value", ylabel = "Count")
    hist!(ax3, vec(Bs), bins = 50)

    return fig
end
fig = sample_and_plot_histograms(3, 3, 3, 10000; lower=-1., upper=1.)