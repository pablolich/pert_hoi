# This file was auto generated from Macaulay2.
# For more information about HomotopyContinuation.jl visit:
#    https://www.JuliaHomotopyContinuation.org
using HomotopyContinuation
using Plots, ImplicitPlots

@var A11, A12, A21, A22, B111, B112, B122, B211, B212, B222, r1, r2

h1 = include("J1.jl")
h2 = include("J2.jl")

@var a

A = [A11, A12, A21, A22]
B = [B111, B112, B122, B211, B212, B222]

g1 = subs(h1, A => (1-a) .* A, B => a .* B)
g2 = subs(h2, A => (1-a) .* A, B => a .* B)

f = System([g1; g2]; variables = [r1;r2], parameters = [a; A; B])

parameter = [2; 4/5; 1/2; 1/8; 4/3; 1; 10/7; 1/3; 4/3; 9/4]
f1 = f([r1;r2], [a; parameter])

l = (-10, 10)

anim = @animate for ai in 0.1:0.1:0.9
    f2 = subs.(f1, a => ai)
    implicit_plot(f2[1]; linecolor = :indianred, linewidth = 2, label = "zero boundary", 
        title = "α = $(ai)", 
        xlabel = "r1", 
        ylabel = "r2",
        xlims = l,
        ylims = l)
    implicit_plot!(f2[2]; linecolor = :steelblue, linewidth = 2, label = "bifurcation",
        xlims = l,
        ylims = l)
end

gif(anim, "anim_fps15.gif", fps = 15)