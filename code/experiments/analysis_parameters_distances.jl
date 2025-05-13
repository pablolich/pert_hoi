using Glob, Serialization

using CSV, DataFrames, Statistics, CategoricalArrays, StatsBase, CairoMakie


function stack_parameter_seed_order(pars::Tuple, seed_i::Int)
    param_rows = Vector{NamedTuple{(:seed_i, :order, :value), Tuple{Int, Int, Float64}}}()

    for (i, param_group) in enumerate(pars)
        order = i - 1
        values = param_group[:]
        append!(param_rows, [(seed_i=seed_i, order=order, value=val) for val in values])
    end

    return DataFrame(param_rows)
end


function load_all_parameters_stack(folder::String, n::Int)
    files = Glob.glob("parameters_n$(n)_seed*.bin", folder)
    all_dfs = DataFrame[]

    for file in files
        seed_i, pars = deserialize(open(file))
        push!(all_dfs, stack_parameter_seed_order(pars, seed_i))
    end

    return vcat(all_dfs...)
end

function histogram_by_order_plot(df::DataFrame)
    fig = Figure(resolution = (1200, 400))
    orders = sort(unique(df.order))
    robustness_classes = ["INCREASE", "DECREASE"]
    colors = (; INCREASE = :blue, DECREASE = :red)

    for (i, order_val) in enumerate(orders)
        ax = Axis(fig[1, i], title = "Order = $order_val", xlabel = "value", ylabel = "count")
        #Makie.xlims!(ax, -2, 2)

        for class in robustness_classes
            subset = df[(df.order .== order_val) .&& (df.robustness .== class), :]
            if nrow(subset) > 0
                hist!(ax, subset.value;
                    bins = 50,
                    normalization = :none,
                    strokewidth = 0.5,
                    color = (colors[Symbol(class)], 0.4),  # transparent fill
                    label = class)
            end
        end

        axislegend(ax, position = :rt)
    end

    return fig
end

df = CSV.read("../../data/boundary_distances/n_3.csv", DataFrame)

rename!(df, [:seed_i, :n, :d, :pert_size, :pert_i, :alpha_i, :δₚ, :δₚₗᵢₙ, :δₓ, :δₓₗᵢₙ, :flag])

# Group by specified columns and compute mean and variance of δₚ
agg_df = combine(DataFrames.groupby(df, [:seed_i, :n, :d, :pert_size, :alpha_i]),
    :δₚ => mean => :δₚ_mean,
    :δₚ => var => :δₚ_var)

# Unstack δₚ_mean by alpha_i
mean_df = unstack(agg_df, [:seed_i, :n, :d, :pert_size], :alpha_i, :δₚ_mean)
rename!(mean_df, names(mean_df, r"^0.1$")[1] => :mean_0_1,
                 names(mean_df, r"^0.9$")[1] => :mean_0_9)

# Unstack δₚ_var by alpha_i
var_df = unstack(agg_df, [:seed_i, :n, :d, :pert_size], :alpha_i, :δₚ_var)
rename!(var_df, names(var_df, r"^0.1$")[1] => :var_0_1,
                 names(var_df, r"^0.9$")[1] => :var_0_9)

# Join and compute deltas
joined = innerjoin(mean_df, var_df, on=[:seed_i, :n, :d, :pert_size])
joined.dmean = joined.mean_0_9 .- joined.mean_0_1
joined.dvar = joined.var_0_9 .- joined.var_0_1

# Select result
result_df = select(joined, :seed_i, :n, :d, :pert_size, :dmean, :dvar)

result_df = transform(joined, [:dmean, :dvar] => ByRow((dm, dv) -> 
    ismissing(dm) || ismissing(dv) ? "UNCLEAR" :
    dm > 0 && dv > 0 ? "INCREASE" :
    dm < 0 && dv < 0 ? "DECREASE" : "UNCLEAR") => :robustness)

result_df

param_df = load_all_parameters_stack("../../data/parameter_sets", 3)

# Join with result_df to bring in robustness label
annotated_params = leftjoin(param_df, result_df, on=:seed_i)

first(annotated_params, 10)  # preview
using StatsPlots, DataFrames
dropmissing!(annotated_params)

# Filter out UNCLEAR and missing values
filtered_df = filter(row -> row.robustness != "UNCLEAR", annotated_params)
dropmissing!(filtered_df)
filtered_df.robustness = categorical(filtered_df.robustness)

fig = histogram_by_order_plot(filtered_df)
save("../../figures/parameter_density_by_robustness_makie.pdf", fig)
