using CSV, DataFrames
using CairoMakie

function plot_spaghetti_by_seed(pairwise_csv;
                                outdir::AbstractString="figures/pairwise_spaghetti",
                                ncols::Int=5,
                                markersize::Real=6,
                                log_eps::Real=1e-12,
                                color_rr=:black,
                                color_rc=:crimson)

    df = CSV.read(pairwise_csv, DataFrame)
    seeds = sort(unique(df.seed))
    mkpath(outdir)

    for seed in seeds
        sub_seed = df[df.seed .== seed, :]
        pert_list = sort(unique(sub_seed.pert_i))
        npert = length(pert_list)
        ncols_eff = max(1, min(ncols, npert))
        nrows = cld(npert, ncols_eff)

        f = Figure(resolution = (ncols_eff * 360, nrows * 260))
        Label(f[0, 1:ncols_eff], "Distances to x*₍crit₎ vs α — seed $(seed)",
              fontsize=18, tellwidth=false)

        axs = Axis[]
        for (k, pert) in enumerate(pert_list)
            r = div(k-1, ncols_eff) + 1
            c = mod(k-1, ncols_eff) + 1
            ax = Axis(f[r, c], title = "pert = $(pert)", yscale = log10)
            push!(axs, ax)

            gpert = sub_seed[sub_seed.pert_i .== pert, :]

            # --- filter rows with missings in plotting-critical columns ---
            keep = .!ismissing.(gpert.alpha) .&
                   .!ismissing.(gpert.dist) .&
                   .!ismissing.(gpert.sol_j_type)
            gpert = gpert[keep, :]

            if nrow(gpert) == 0
                text!(ax, 0.5, 0.5, text = "no pairs",
                      align = (:center, :center), space = :relative)
                continue
            end

            # log-scale safe y values
            yvals = max.(gpert.dist, log_eps)

            # masks (now pure Bool; no missings after filter)
            mask_real    = gpert.sol_j_type .== "real"
            mask_complex = gpert.sol_j_type .== "complex"

            h_objs = Any[]
            labels = String[]

            if any(mask_real)
                h_rr = scatter!(ax, gpert.alpha[mask_real], yvals[mask_real];
                                marker = :circle, color = color_rr,
                                markersize = markersize)
                push!(h_objs, h_rr); push!(labels, "real–real")
            end
            if any(mask_complex)
                h_rc = scatter!(ax, gpert.alpha[mask_complex], yvals[mask_complex];
                                marker = :circle, color = color_rc,
                                markersize = markersize)
                push!(h_objs, h_rc); push!(labels, "real–complex")
            end

            # y-lims: start fixed at 1e-4, ensure upper > lower
            ymin = 1e-4
            ymax = maximum(yvals)
            ymax_plot = max(ymax * 1.3, ymin * 1.1)
            ylims!(ax, ymin, ymax_plot)

            # Legend only on the first panel
            if k == 1 && !isempty(h_objs)
                axislegend(ax, h_objs, labels; position = :rt, framevisible = false)
            end
        end

        # label outer axes only
        for c in 1:ncols_eff
            axs[(nrows - 1) * ncols_eff + c].xlabel = "α"
        end
        for r in 1:nrows
            axs[(r - 1) * ncols_eff + 1].ylabel = "distance to x*₍crit₎"
        end

        if length(axs) > 1
            linkxaxes!(axs...)
        end

        outpath = joinpath(outdir, "pairwise_spaghetti_seed_$(seed).png")
        CairoMakie.save(outpath, f)
        @info "Saved" outpath
    end

    nothing
end


ns = [2, 3]
for n in ns
    csv = joinpath("simulation_results", "solution_pairwise_distances_n_$(n).csv")
    isfile(csv) || (@warn "Missing CSV, skipping" csv; continue)
    outdir = joinpath("figures/pairwise_spaghetti", "n_$(n)")
    mkpath(outdir)
    plot_spaghetti_by_seed(csv; outdir=outdir, ncols=5, markersize=6)
end
