using HomotopyContinuation

@var x a y b




F = System([x^2 - a], [x], [a])

ct = Tracker(

ParameterHomotopy(F, [1], [2]),

options = TrackerOptions(max_step_size = 0.015625),

)

Xs = Vector{ComplexF64}[]

for (x, t) in iterator(ct, [-1.0], 1.0, 0.0)

push!(Xs, x)

end

