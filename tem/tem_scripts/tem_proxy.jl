# Here I will implement my TEM-proxy Mathematica scripts into Julia

using DelimitedFiles, CairoMakie, GLMakie, LaTeXStrings, Trapz

# Extract mode structures from zprofiles  

zp1 = readdlm("zprofileions_0001.dat", skipstart=0)

ϕ1 = zp1[:, 3];
ϕ1 = convert(Array{Float64}, ϕ1);

nz = length(ϕ1);

ϕ1mx = maximum(ϕ1);
ϕ1nrm = ϕ1./ϕ1mx

ϕ1 = ϕ1nrm

# Extract geometry from gist: B, gyy, κ

gist = readdlm("gist.dat", skipstart=0);

B = gist[:,4];
B = convert(Array{Float64}, B);
Bmax = maximum(B);
Bmin = minimum(B);

gyy = gist[:,3];
gyy = convert(Array{Float64}, gyy);

κ = gist[:,4];
κ = convert(Array{Float64}, κ);

# Pitch angle λ

λ = 1/Bmax:0.001:1/Bmin

well = heaviside(1/λ - B)

l = range(1, nz, length=nz)
τ = trapz(well/sqrt(1-λ*B)

