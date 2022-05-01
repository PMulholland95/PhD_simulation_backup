# Here I will implement my TEM-proxy Mathematica scripts into Julia

using DelimitedFiles, Trapz, Kinetic, LaTeXStrings, SpecialFunctions

# Extract geometry from gist: B, gyy, κ
# Safety factors: d3d, ncsx, w7xsc, w7xhm, w7xlm, kjm, dkh, dkm, dks

function omegaTemProxyPhiGene(gistread::String,
		       	      scanlog::String,
			      zprofs::String,
		       	      savefolder::String,
			      vn::Int64,
			      sf::Int64)

	n = [11,17]
	zn= [44,68]

	gist = readdlm(gistread, skipstart=13);

	B = gist[:,4];
	B = convert(Array{Float64}, B);
	Bmax = maximum(B);
	Bmin = minimum(B);

	writedlm(string(savefolder,"B_z.dat"),B)

	gyy = gist[:,3];
	gyy = convert(Array{Float64}, gyy);

	writedlm(string(savefolder,"gyy_z.dat"),gyy)

	kapn = [5,6]	

	κ = gist[:,kapn[vn]];
	κ = convert(Array{Float64}, κ);

	writedlm(string(savefolder,"kappa_z.dat"),κ)

	mydpdx = 0;
	s = 0.5;

	# Importing GENE mode structure

	s1 = ["zprofileions_000","zprofilei_000"]
	s3 = ".dat"

	v = []

	for i in 1:9
		push!(v, string(s1[vn], i, s3))
	end

	s2 = ["zprofileions_00","zprofilei_00"] 

	for i in 10:zn[vn]
		push!(v, string(s2[vn], i, s3))
	end

	sc = []
	ϕall0 = []

	for i in 1:zn[vn]
		push!(sc, readdlm(string(zprofs,v[i]),skipstart=7))
		push!(ϕall0, sc[i][:,3])
	end

	ϕall = []

	for i in 1:zn[vn]
		push!(ϕall, ϕall0[i] ./ maximum(ϕall0[i]))
	end
	
	pmatempty = Array{Vector{Float64}}(undef, n[vn], 4)

	for i in eachindex(pmatempty)
	pmatempty[i] = ϕall[i]
	end 

	pmat = pmatempty

	z = sc[1][:,1]

	# Mode Frequency Proxy

	λ = range(1/Bmax, 1/Bmin, length=1000);
	l = 1:1:length(B);
	m = 1:1:zn[vn];

	kyall = [[0.6,0.8,1.0,1.2], [0.8,1.0,1.2,1.4]]
	gradall = [[0.3,0.6,0.8,1.0,1.2,1.5,2.0,3.0,4.0,5.0,6.0], [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0]]

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

		kn = indexin(kyρ, kyall[vn])
		gn = indexin(gradn, gradall[vn])

		(1 ./BounceTime(λ)).*real(trapz((l), (Well(λ).* pmat[gn[1],kn[1]])./sqrt.(Complex.(1 .- λ .* B)))) 
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

		kn = indexin(kyρ, kyall[vn])
		gn = indexin(gradn, gradall[vn])

		(trapz((l), (2 .- Γ0(kyρ)) .* (pmat[gn[1],kn[1]] .^2) .* (1 ./B)))
	end

	function a2(kyρ, gradn)
		
		0.5 .* (trapz((λ), BIL(kyρ, gradn)))
	end

	function b1(kyρ, gradn, ηᵢ)

		kn = indexin(kyρ, kyall[vn])
		gn = indexin(gradn, gradall[vn])

		(trapz((l), (Γ0(kyρ) .- ηᵢ * b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* (pmat[gn[1],kn[1]] .^2) .* (1 ./B))) 
	end

	function b2(kyρ, gradn)

		kn = indexin(kyρ, kyall[vn])
		gn = indexin(gradn, gradall[vn])

		(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* κ .* (pmat[gn[1],kn[1]] .^2) .* (1 ./B))) 
	end

	function b3(kyρ, gradn)

		(1/gradn)*(-3/2) * (trapz((λ), BIL(kyρ, gradn) .* GBL))
	end

	function c1(kyρ, gradn, ηᵢ)
		
		kn = indexin(kyρ, kyall[vn])
		gn = indexin(gradn, gradall[vn])

		(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ)) + ηᵢ .* (2 .* ((b(kyρ) .-1) .^2) .* Γ0(kyρ) .+ b(kyρ) .* (3 .- 2 .* b(kyρ)) .* Γ1(kyρ) )) .* κ .* (pmat[gn[1],kn[1]] .^2) .* (1 ./B) )) 
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

	gsc = readdlm(scanlog, skipstart=1);

	# Safety factors: d3d, ncsx, w7xsc, w7xhm, w7xlm, kjm, dkh, dkm, dks
	q0 = [2.5655027, 1.8588439, 1.1185532, 1.101862399, 1.1266741, 1.10964797, 1.1086111079, 1.092305657, 1.0852873343]  

	g1 = gsc[1:n[vn],8]
	g2 = gsc[n[vn]+1:2n[vn],8]
	g3 = gsc[2n[vn]+1:3n[vn],8]
	g4 = gsc[3n[vn]+1:4n[vn],8]

	gradn = gsc[1:n[vn],3]
	gradn = convert(Array{Float64}, gradn)

	z1 = zeros(n[vn],1)
	z2 = zeros(n[vn],1)
	z3 = zeros(n[vn],1)
	z4 = zeros(n[vn],1)

	k1 = [0.6,0.8]
	k2 = [0.8,1.0]
	k3 = [1.0,1.2]
	k4 = [1.2,1.4]

	ky1 = fill!(z1,k1[vn])
	ky2 = fill!(z2,k2[vn])
	ky3 = fill!(z3,k3[vn])
	ky4 = fill!(z4,k4[vn])

	etai = zeros(n[vn])
	etae = zeros(n[vn])

	qp1 = map(x->quadplus(ky1[x],gradn[x],etai[x],etae[x]), 1:n[vn])
	qp2 = map(x->quadplus(ky2[x],gradn[x],etai[x],etae[x]), 1:n[vn])
	qp3 = map(x->quadplus(ky3[x],gradn[x],etai[x],etae[x]), 1:n[vn])
	qp4 = map(x->quadplus(ky4[x],gradn[x],etai[x],etae[x]), 1:n[vn])

	qpr1 = real(qp1)
	qpr2 = real(qp2)
	qpr3 = real(qp3)
	qpr4 = real(qp4)

	ban1 = vec(-(q0[sf] .* ky1 .* gradn) .* qpr1) 
	ban2 = vec(-(q0[sf] .* ky2 .* gradn) .* qpr2) 
	ban3 = vec(-(q0[sf] .* ky3 .* gradn) .* qpr3) 
	ban4 = vec(-(q0[sf] .* ky4 .* gradn) .* qpr4) 

	qpr = "qpr"
	ban = "ban"
	g = "g"
	dat = ".dat"
	
	qprar = []
	banar = []
	gar = []

	for i in 1:4
		push!(qprar, string(savefolder,qpr,i,dat))
		push!(banar, string(savefolder,ban,i,dat))
		push!(gar, string(savefolder,g,i,dat))
	end
	
	writedlm(qprar[1],qpr1)
	writedlm(qprar[2],qpr2)
	writedlm(qprar[3],qpr3)
	writedlm(qprar[4],qpr4)

	writedlm(banar[1],ban1)
	writedlm(banar[2],ban2)
	writedlm(banar[3],ban3)
	writedlm(banar[4],ban4)

	writedlm(gar[1],g1)
	writedlm(gar[2],g2)
	writedlm(gar[3],g3)
	writedlm(gar[4],g4)

end
