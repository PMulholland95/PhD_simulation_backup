nz = 1:1:64

function Well(λ::Float64)

	map(x->heaviside(1/λ - B[x]),nz)
end

# τ = trapz((l), well/sqrt.(Complex.(1 .- λ.*B)))

function BounceTime(λ::Float64)

	trapz((l), Well(λ)./sqrt.(Complex.(1 .- λ*B)))
end


