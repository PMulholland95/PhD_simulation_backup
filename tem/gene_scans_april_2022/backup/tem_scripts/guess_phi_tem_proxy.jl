# Here I will implement my TEM-proxy Mathematica scripts into Julia

using DelimitedFiles, Trapz, Kinetic, LaTeXStrings, SpecialFunctions

# Extract geometry from gist: B, gyy, κ

function OmegaTemProxyPhiGuess(gistread::String,
		       	       scanlog::String,
		       	       savefolder::String)
		       

	gist = readdlm(gistread, skipstart=0);

	B = gist[:,4];
	B = convert(Array{Float64}, B);
	Bmax = maximum(B);
	Bmin = minimum(B);

	gyy = gist[:,3];
	gyy = convert(Array{Float64}, gyy);

	κ = gist[:,5];
	κ = convert(Array{Float64}, κ);

	mydpdx = 0;
	q = 1;
	s = 0.5;

	# Predicting the mode structure

	bk = B + κ;
	inbk = 1 ./bk;
	inbkmx = maximum(inbk);

	ϕGuessNorm = inbk/inbkmx;

	ϕ = ϕGuessNorm;

	# Mode Frequency Proxy

	# λ = 1/Bmax:0.001:1/Bmin
	λ = range(1/Bmax, 1/Bmin, length=1000);
	l = 1:1:length(B);

	function Well(λ::Float64)

		map(x->heaviside(1/λ - B[x]),l)
	end

	well = map(x->Well(x), λ)

	# τ = trapz((l), well/sqrt.(Complex.(1 .- λ.*B)))

	function BounceTime(λ::Float64)

		real(trapz((l), Well(λ)./sqrt.(Complex.(1 .- λ*B))))
	end

	τ = map(x->BounceTime(x), λ);

	# τreal = real(τ);
	# τimag = imag(τ);

	function PhiBounce(λ::Float64)

		(1 ./BounceTime(λ)).*real(trapz((l), (Well(λ).*ϕ)./sqrt.(Complex.(1 .- λ*B)))) 
	end

	ϕBounceList = map(x->PhiBounce(x), λ);
	ϕBL = ϕBounceList

	BounceIntList = (ϕBL.^2).*τ 
	BIL = BounceIntList
	replace!(BIL, NaN => 0)


	# Procedure for writing ω_d (bounce part)

	function Gbounce(λ::Float64)

		(1 ./BounceTime(λ)).*real(trapz((l), (Well(λ).*κ).*(1 .- (λ*B./2))./sqrt.(Complex.(1 .- λ*B))))
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

	function a1(kyρ::Float64)

		(trapz((l), (2 .- Γ0(kyρ)).*(ϕ.^2).*(1 ./B)))
	end

	a2 = 0.5 .* (trapz((λ), BIL))

	function b1(kyρ, ηᵢ)

		(trapz((l), (Γ0(kyρ) .- ηᵢ * b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* (ϕ .^2) .* (1 ./B))) 
	end

	function b2(kyρ, gradn)

		(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ))) .* κ .* (ϕ .^2) .* (1 ./B))) 
	end

	function b3(kyρ, gradn)

		(1/gradn)*(-3/2) * (trapz((λ), BIL .* GBL))
	end

	function c1(kyρ, gradn, ηᵢ)
		
		(1/gradn) * (trapz((l), (2 * Γ0(kyρ) .- b(kyρ) .* (Γ0(kyρ) - Γ1(kyρ)) + ηᵢ .* (2 .* ((b(kyρ) .-1) .^2) .* Γ0(kyρ) .+ b(kyρ) .* (3 .- 2 .* b(kyρ)) .* Γ1(kyρ) )) .* κ .* (ϕ .^2) .* (1 ./B) )) 
	end

	function c2(kyρ, gradn, ηₑ)

		b3(kyρ, gradn) .* (1 .+ ηₑ)
	end

	# Grouping together into quadratic terms

	function aquad(kyρ)

		a1(kyρ) .- a2
	end

	function bquad(kyρ, gradn, ηᵢ)

		b1(kyρ, ηᵢ) .+ b2(kyρ, gradn) .- a2 .+ b3(kyρ, gradn)
	end

	function cquad(kyρ, gradn, ηᵢ, ηₑ)

		c1(kyρ, gradn, ηᵢ) - c2(kyρ, gradn, ηₑ)
	end

	function quadplus(kyρ, gradn, ηᵢ, ηₑ)

		(bquad(kyρ, gradn, ηᵢ) .+ sqrt.(Complex.((bquad(kyρ, gradn, ηᵢ)).^2 + 4 .* aquad(kyρ) .* cquad(kyρ, gradn, ηᵢ, ηₑ)))) ./ (2 .* aquad(kyρ))

	end

	# Now, import relevant data from scan.log

	gsc = readdlm(scanlog);

	n = 11
	q0 = 2.5655027

	g1 = gsc[1:n,8]
	g2 = gsc[n+1:2n,8]
	g3 = gsc[2n+1:3n,8]
	g4 = gsc[3n+1:4n,8]

	gradn = gsc[1:n,3]
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
	writedlm(qprar[4],qpr3)

	writedlm(banar[1],ban1)
	writedlm(banar[2],ban2)
	writedlm(banar[3],ban3)
	writedlm(banar[4],ban3)

	writedlm(gar[1],g1)
	writedlm(gar[2],g2)
	writedlm(gar[3],g3)
	writedlm(gar[4],g3)

end
