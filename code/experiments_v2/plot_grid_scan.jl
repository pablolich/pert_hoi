# plot_grid_scan.jl
#
# Load a serialized grid scan artifact created by compute_grid_scan.jl and
# produce two PNGs:
#   (A) multi-panel heatmaps of min real component with boundary overlay
#   (B) multi-panel heatmaps of min abs imaginary component with boundary overlay

using CairoMakie, Serialization, FilePathsBase

# Flag → color mapping for boundary overlay (tune as you like)
const FLAG_COLORS = Dict(
    "baseline"     => "#377EB8",
    "negative"     => "#377EB8",
    "complex"      => "#984EA3",
    "boundary"     => "#4DAF4A",
    "nonconverged" => "#E41A1C",
    "unknown"      => "black",
)

# Make a small legend for used flags on the first panel in each figure
function _legend_for_flags!(ax, used::Vector{String})
    hs, labs = Any[], String[]
    for k in used
        h = scatter!(ax, [-Inf], [-Inf]; color=FLAG_COLORS[k], markersize=8)
        push!(hs, h); push!(labs, k)
    end
    if !isempty(hs)
        axislegend(ax, hs, labs; position=:rt, framevisible=false)
    end
end

function plot_grid_scan(inpath::AbstractString; outdir::AbstractString="grid_scan_plots",
                        ncols::Int=3, fontsize::Int=14)

    isdir(outdir) || mkpath(outdir)

    data = open(inpath, "r") do io
        deserialize(io)
    end

    n            = data[:n]
    seed         = data[:seed]
    n_side       = data[:n_side]
    pert_size    = data[:pert_size]
    alphas       = data[:alphas]
    x_axis       = data[:x_axis]
    y_axis       = data[:y_axis]
    zreal_slices = data[:zreal_slices]
    zimag_slices = data[:zimag_slices]
    boundary_pts = data[:boundary_pts]
    boundary_flags = data[:boundary_flags]

    nα    = length(alphas)
    ncols = min(ncols, max(1, nα))
    nrows = cld(nα, ncols)

    # -------------------------------
    # (A) min real component heatmaps
    # -------------------------------
    figA = Figure(resolution=(400*ncols, 360*nrows), fontsize=fontsize)
    axsA = Axis[]

    # symmetric clims for real (center color at 0)
    # compute global symmetric clims to make panels comparable
    zmin = minimum(map(minimum, zreal_slices))
    zmax = maximum(map(maximum, zreal_slices))
    clim_abs = max(abs(zmin), abs(zmax))
    clims_real = (-clim_abs, clim_abs)

    for k in 1:nα
        r = div(k-1, ncols) + 1
        c = mod(k-1, ncols) + 1
        ax = Axis(figA[r, c], title = "α = $(round(alphas[k], digits=3))")
        push!(axsA, ax)

        hm = heatmap!(ax, x_axis, y_axis, zreal_slices[k];
                      colormap = cgrad(:RdBu, rev=true), colorrange = clims_real)
        Colorbar(figA[r, c+1], hm, label="min(real component)"; height=Relative(0.8))

        # overlay boundary points
        pts   = boundary_pts[k]
        flags = boundary_flags[k]
        used  = String[]
        if size(pts, 1) > 0
            for (j, lab) in enumerate(flags)
                color = get(FLAG_COLORS, lab, "black")
                scatter!(ax, [pts[j,1]], [pts[j,2]];
                         color=color, markersize=5)
                push!(used, lab)
            end
        end

        # axes labels on outer only
        if r == nrows; ax.xlabel = "r₁"; end
        if c == 1;     ax.ylabel = "r₂"; end

        if k == 1 && !isempty(used)
            _legend_for_flags!(ax, unique(used))
        end
    end

    outA = joinpath(outdir, "heatmap_minreal_n_$(n)_seed_$(seed).png")
    CairoMakie.save(outA, figA)

    # # -----------------------------------------
    # # (B) min abs imaginary component heatmaps
    # # -----------------------------------------
    # figB = Figure(resolution=(400*ncols, 360*nrows), fontsize=fontsize)
    # axsB = Axis[]

    # # sequential clims for abs(imag): nonnegative
    # zminI = minimum(map(minimum, zimag_slices))
    # zmaxI = maximum(map(maximum, zimag_slices))
    # clims_imag = (zminI, zmaxI)

    # for k in 1:nα
    #     r = div(k-1, ncols) + 1
    #     c = mod(k-1, ncols) + 1
    #     ax = Axis(figB[r, c], title = "α = $(round(alphas[k], digits=3))")
    #     push!(axsB, ax)

    #     hm = heatmap!(ax, x_axis, y_axis, zimag_slices[k];
    #                   colormap = :viridis, colorrange = clims_imag)
    #     Colorbar(figB[r, c+1], hm, label="min |imag component|"; height=Relative(0.8))

    #     # overlay boundary points
    #     pts   = boundary_pts[k]
    #     flags = boundary_flags[k]
    #     used  = String[]
    #     if size(pts, 1) > 0
    #         for (j, lab) in enumerate(flags)
    #             color = get(FLAG_COLORS, lab, "black")
    #             scatter!(ax, [pts[j,1]], [pts[j,2]];
    #                      color=color, markersize=5)
    #             push!(used, lab)
    #         end
    #     end

    #     if r == nrows; ax.xlabel = "r₁"; end
    #     if c == 1;     ax.ylabel = "r₂"; end

    #     if k == 1 && !isempty(used)
    #         _legend_for_flags!(ax, unique(used))
    #     end
    # end

    # outB = joinpath(outdir, "heatmap_minimag_n_$(n)_seed_$(seed).png")
    #CairoMakie.save(outB, figB)

    println("Saved:")
    println("  → ", outA)
    #println("  → ", outB)
end

# ---------
# Example(s)
# ---------
plot_grid_scan("grid_scans/grid_scan_n_2_seed_3.bin"; outdir="grid_scan_plots", ncols=3)
