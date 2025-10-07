using CairoMakie, CSV, DataFrames
using FilePathsBase: mkpath, joinpath

n_values = [2]

flag_colors = Dict(
    "baseline"     => "#377EB8",
    "negative"     => "#377EB8",
    "complex"      => "#984EA3",
    "boundary"     => "#4DAF4A",
    "nonconverged" => "#E41A1C",
    "unknown"      => "black"
)

flag_shapes = Dict(
    "baseline"     => :hline,
    "negative"     => :hline,
    "complex"      => :star4,
    "boundary"     => :cross,
    "nonconverged" => :xcross,
    "unknown"      => :star8
)

outdir = "plots_many_alphas"
mkpath(outdir)

for n in n_values
    println("Processing n = $n ...")

    infile = "simulation_results/boundary_distances_n_$(n)_many_alphas_many_perts.csv"
    df = CSV.read(infile, DataFrame)
    rename!(df, Symbol("δ_p") => :dp, Symbol("δ_x") => :dx)

    seeds    = unique(df.seed)
    perturbs = unique(df.pert_i)

    for seed in seeds
        subdf = df[df.seed .== seed, :]

        fig = Figure(resolution = (900, 300))  # wider for legend
        grid = fig[1, 1:2] = GridLayout()

        ax1 = Axis(grid[1, 1];
            xlabel = "α",
            ylabel = "δₚ",
            title  = "n = $n — Parameter Distance",
            yscale = log10)

        ax2 = Axis(grid[1, 2];
            xlabel = "α",
            ylabel = "δₓ/δₚ",
            title  = "n = $n — Equilibrium Distance",
            yscale = log10)

        # Store one scatter! handle per flag (for legend)
        legend_handles = Dict{String, Any}()

        for pert in perturbs
            d = subdf[subdf.pert_i .== pert, :]

            lines!(ax1, d.alpha, d.dp, color = :gray, linewidth = 1)
            lines!(ax2, d.alpha, d.dx, color = :gray, linewidth = 1)

            # Per-point style
            flags = string.(d.flag)
            colors = [get(flag_colors, f, "black") for f in flags]
            shapes = [get(flag_shapes, f, :circle) for f in flags]

            # Scatter each point individually so we can track one representative per flag
            for (α, δp, δx, c, m, f) in zip(d.alpha, d.dp, d.dx, colors, shapes, flags)
                h1 = scatter!(ax1, [α], [δp], color=c, marker=m, markersize=7)
                scatter!(ax2, [α], [δx], color=c, marker=m, markersize=7)

                # Store the first handle per flag for legend
                if !haskey(legend_handles, f)
                    legend_handles[f] = h1
                end
            end
        end

        # Filter out "baseline" from legend (same as "negative")
        legend_keys = filter(k -> k != "baseline", keys(legend_handles))
        sorted_keys = sort(collect(legend_keys))
        handles = [legend_handles[k] for k in sorted_keys]
        labels = sorted_keys


        # Add legend outside to the right
        fig[1, 3] = Legend(fig, handles, labels, "Flag Type"; tellheight=false)

        outfile = joinpath(outdir, "n_$(n)_seed_$(seed).png")
        CairoMakie.save(outfile, fig)
    end

    println("✅ Finished plots for n = $n")
end

println("🎉 All plots saved to: $outdir")
