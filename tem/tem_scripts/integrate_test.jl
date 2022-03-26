using Interpolations, CairoMakie, Trapz

# Example data to interpolate
x = LinRange(0,π,50)
z(x) = sin(x)
y = z.(x)

M = Array(y)

I = trapz((x),M)


