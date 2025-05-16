# This file was auto generated from Macaulay2.
# For more information about HomotopyContinuation.jl visit:
#    https://www.JuliaHomotopyContinuation.org
using HomotopyContinuation
using Plots, ImplicitPlots

@var r1 r2 a 

h1 = include("J1.jl")
h2 = include("J2.jl")
f = System([h1; h2]; variables = [r1;r2], parameters = [a])

f1 = f([r1;r2], [0.1])
f2 = f([r1;r2], [0.9])

l = (-10, 10)
implicit_plot(f1[1]; linecolor = :indianred, linewidth = 2, label = "zero boundary", 
    title = "α = 0.1", 
    xlabel = "r1", 
    ylabel = "r2",
    xlims = l,
    ylims = l)
implicit_plot!(f1[2]; linecolor = :steelblue, linewidth = 2, label = "bifurcation",
    xlims = l,
    ylims = l)
scatter!([1],[1], markercolor = :black, markersize = 3, label = "(1,1)")


implicit_plot(f2[1]; linecolor = :indianred, linewidth = 2, label = "zero boundary", title = "α = 0.9", 
    xlabel = "r1", 
    ylabel = "r2",
    xlims = l,
    ylims = l)
implicit_plot!(f2[2]; linecolor = :steelblue, linewidth = 2, label = "bifurcation",
    xlims = l,
    ylims = l)
scatter!([1],[1], markercolor = :black, markersize = 3, label = "(1,1)")