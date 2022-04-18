# Here I will implement my TEM-proxy Mathematica scripts into Julia

using DelimitedFiles, Trapz, Kinetic, LaTeXStrings, SpecialFunctions

# Extract geometry from gist: B, gyy, κ

gist = readdlm("start_data/gist_d3d_test.dat", skipstart=0);

B = gist[:,4];
B = convert(Array{Float64}, B);
Bmax = maximum(B);
Bmin = minimum(B);

gyy = gist[:,3];
gyy = convert(Array{Float64}, gyy);

κ = gist[:,5];
κ = convert(Array{Float64}, κ);

mydpdx = 0;
s = 0.5;

# Importing GENE mode structure

s1 = "zprofileions_000"
s3 = ".dat"

v = []

for i in 1:9
	push!(v, string(s1, i, s3))
end

s2 = "zprofileions_00" 

for i in 10:44
	push!(v, string(s2, i, s3))
end

sc = []
ϕall0 = []

for i in 1:44
	push!(sc, readdlm(string("start_data/d3d_phi_profiles/",v[i]),skipstart=7))
	push!(ϕall0, sc[i][:,3])
end

ϕall = []

for i in 1:44
	push!(ϕall, ϕall0[i] ./ maximum(ϕall0[i]))
end

pmat = [[ϕall[1],ϕall[12],ϕall[23],ϕall[34]] [ϕall[2],ϕall[13],ϕall[24],ϕall[35]] [ϕall[3],ϕall[14],ϕall[25],ϕall[36]] [ϕall[4],ϕall[15],ϕall[26],ϕall[37]] [ϕall[5],ϕall[16],ϕall[27],ϕall[38]] [ϕall[6],ϕall[17],ϕall[28],ϕall[39]] [ϕall[7],ϕall[18],ϕall[29],ϕall[40]] [ϕall[8],ϕall[19],ϕall[30],ϕall[41]] [ϕall[9],ϕall[20],ϕall[31],ϕall[42]] [ϕall[10],ϕall[21],ϕall[32],ϕall[43]] [ϕall[11],ϕall[22],ϕall[33],ϕall[44]]]

z = sc[1][:,1]

# Mode Frequency Proxy

# λ = 1/Bmax:0.001:1/Bmin
λ = range(1/Bmax, 1/Bmin, length=1000);
l = 1:1:length(B);
m = 1:1:44;
kyall = [0.6,0.8,1.0,1.2]
gradall = [0.3,0.6,0.8,1.0,1.2,1.5,2.0,3.0,4.0,5.0,6.0]

function Well(λ::Float64)

	map(x->heaviside(1 ./λ .- B[x]),l)
end

well = map(x->Well(x), λ)

# τ = trapz((l), well/sqrt.(Complex.(1 .- λ.*B)))

function BounceTime(λ::Float64)

	real(trapz((l), Well(λ)./sqrt.(Complex.(1 .- λ .* B))))
end

τ = map(x->BounceTime(x), λ);
replace!(τ, NaN => 0)

# τreal = real(τ);
# τimag = imag(τ);

function PhiBounce(λ::Float64, kyρ, gradn)

	kn = indexin(kyρ, kyall)
	gn = indexin(gradn, gradall)
	kg = kn .* gn

	(1 ./BounceTime(λ)).*real(trapz((l), (Well(λ).*pmat[kn[1],gn[1]])./sqrt.(Complex.(1 .- λ .* B)))) 
end

function ϕBL(kyρ, gradn)

	ϕBounceList = map(x->PhiBounce(x, kyρ, gradn), λ)
end

function BIL(kyρ, gradn)

	BounceIntList = (ϕBL(kyρ, gradn) .^2).*τ
	replace!(BounceIntList, NaN => 0)
end

# Procedure for writing ω_d (bounce part)

function Gbounce(λ::Float64)

	(1 ./BounceTime(λ)).*real(trapz((l), (Well(λ).*κ).*(1 .- (λ .* B ./2))./sqrt.(Complex.(1 .- λ .* B))))
end

GBounceList = map(x->Gbounce(x), λ);
GBL = GBounceList
replace!(GBL, NaN => 0)

# Postive quadratic soluiton for VP ω-proxy

function b(kyρ::Float64)

	gyy .*(kyρ^2)	
end

function Γ0(kyρ::Float64) 
	
	besseli.(0, b(kyρ)) .* exp.(-b(kyρ))
end

function Γ1(kyρ::Float64) 
	
	besseli.(1, b(kyρ)) .* exp.(-b(kyρ))
end

function a1(kyρ, gradn)

	kn = indexin(kyρ, kyall)
	gn = indexin(gradn, gradall)
	kg = kn .* gn

	(trapz((l), (2 .- Γ0(kyρ)).*(pmat[kn[1],gn[1]] .^2).*(1 ./B)))
end

function a2(kyρ, gradn)
	
	0.5 .* (trapz((λ), BIL(kyρ, gradn)))
end

function b1(kyρ, gradn, ηᵢ)

	kn = indexin(kyρ, kyall)
	gn = indexin(gradn, gradall)
	kg = kn .* gn

	(trapz((l), (Γ0(kyρ) .- ηᵢ * b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* (pmat[kn[1],gn[1]] .^2) .* (1 ./B))) 
end

function b2(kyρ, gradn)

	kn = indexin(kyρ, kyall)
	gn = indexin(gradn, gradall)
	kg = kn .* gn

	(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* κ .* (pmat[kn[1],gn[1]] .^2) .* (1 ./B))) 
end

function b3(kyρ, gradn)

	(1/gradn)*(-3/2) * (trapz((λ), BIL(kyρ, gradn) .* GBL))
end

function c1(kyρ, gradn, ηᵢ)
	
	kn = indexin(kyρ, kyall)
	gn = indexin(gradn, gradall)
	kg = kn .* gn

	(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ)) + ηᵢ .* (2 .* ((b(kyρ) .-1) .^2) .* Γ0(kyρ) .+ b(kyρ) .* (3 .- 2 .* b(kyρ)) .* Γ1(kyρ) )) .* κ .* (pmat[kn[1],gn[1]] .^2) .* (1 ./B) )) 
end

function c2(kyρ, gradn, ηₑ)

	b3(kyρ, gradn) .* (1 .+ ηₑ)
end

# Grouping together into quadratic terms

function aquad(kyρ, gradn)

	a1(kyρ, gradn) .- a2(kyρ, gradn)
end

function bquad(kyρ, gradn, ηᵢ)

	b1(kyρ, gradn, ηᵢ) .+ b2(kyρ, gradn) .- a2(kyρ, gradn) .+ b3(kyρ, gradn)
end

function cquad(kyρ, gradn, ηᵢ, ηₑ)

	c1(kyρ, gradn, ηᵢ) .- c2(kyρ, gradn, ηₑ)
end

function quadplus(kyρ, gradn, ηᵢ, ηₑ)

	(bquad(kyρ, gradn, ηᵢ) .+ sqrt.(Complex.((bquad(kyρ, gradn, ηᵢ)).^2 + 4 .* aquad(kyρ, gradn) .* cquad(kyρ, gradn, ηᵢ, ηₑ)))) ./ (2 .* aquad(kyρ, gradn))

end

# Now, import relevant data from scan.log

dsc = readdlm("start_data/d3d_scan.log");

n = 11
q0 = 2.5655027

d1 = dsc[1:n,8]
d2 = dsc[n+1:2n,8]
d3 = dsc[2n+1:3n,8]
d4 = dsc[3n+1:4n,8]

gradn = dsc[1:n,3]
gradn = convert(Array{Float64}, gradn)

z1 = zeros(n,1)
z2 = zeros(n,1)
z3 = zeros(n,1)
z4 = zeros(n,1)

ky1 = fill!(z1,0.6)
ky2 = fill!(z2,0.8)
ky3 = fill!(z3,1.0)
ky4 = fill!(z4,1.2)

etai = zeros(n)
etae = zeros(n)

qp1 = map(x->quadplus(ky1[x],gradn[x],etai[x],etae[x]), 1:n)
qp2 = map(x->quadplus(ky2[x],gradn[x],etai[x],etae[x]), 1:n)
qp3 = map(x->quadplus(ky3[x],gradn[x],etai[x],etae[x]), 1:n)
qp4 = map(x->quadplus(ky4[x],gradn[x],etai[x],etae[x]), 1:n)

qpr1 = real(qp1)
qpr2 = real(qp2)
qpr3 = real(qp3)
qpr4 = real(qp4)

ban1 = vec(-(q0 .* ky1 .* gradn) .* qpr1) 
ban2 = vec(-(q0 .* ky2 .* gradn) .* qpr2) 
ban3 = vec(-(q0 .* ky3 .* gradn) .* qpr3) 
ban4 = vec(-(q0 .* ky4 .* gradn) .* qpr4) 

writedlm("saved_data/save_phi_gene/qpr1.dat",qpr1)
writedlm("saved_data/save_phi_gene/qpr2.dat",qpr2)
writedlm("saved_data/save_phi_gene/qpr3.dat",qpr3)
writedlm("saved_data/save_phi_gene/qpr4.dat",qpr4)

writedlm("saved_data/save_phi_gene/gban1.dat",ban1)
writedlm("saved_data/save_phi_gene/gban2.dat",ban2)
writedlm("saved_data/save_phi_gene/gban3.dat",ban3)
writedlm("saved_data/save_phi_gene/gban4.dat",ban4)

writedlm("saved_data/save_phi_gene/gd1.dat",d1)
writedlm("saved_data/save_phi_gene/gd2.dat",d2)
writedlm("saved_data/save_phi_gene/gd3.dat",d3)
writedlm("saved_data/save_phi_gene/gd4.dat",d4)

#=
fig = Figure(resolution=(4000,2000))

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])

lines!(ax1, gradn, d1, color = :blue)
lines!(ax1, gradn, ban1, color = :red)

save("phi_proxy_test.png",fig)

fig
=#