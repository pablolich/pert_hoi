#some code to understand the functionality of evaluate homotopies

using HomotopyContinuation

@var x y a b c

f = [(2 * x^2 + b^2 * y^3 + 2 * a * x * y)^3, (a + c)^4 * x + y^2]
F = System(f, [x, y], [a, b, c])

p = [5.2, -1.3, 9.3]
q = [2.6, 3.3, 2.3]

@var tvar
h = Homotopy(
    f([x, y] => [x, y], [a, b, c] => tvar * p + (1 - tvar) * q),
    [x, y],
    tvar,
)

homotopy = ParameterHomotopy(F, p, q; compile = false)

m, n = size(homotopy)
u = zeros(ComplexF64, m)
x = randn(ComplexF64, n)
t = randn(ComplexF64)

evaluate!(u, homotopy, x, t)