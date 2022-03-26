using CairoMakie
using LaTeXStrings


function plotModeSpectrum(model::PlasmaTurbulenceSaturationModel.FluidModel;
		resolution=(1600,900),
		fontsize=25,
		title="Mode spectrum",
		)

	kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
 	ky = PlasmaTurbulenceSaturationModel.ky(model)

	evals = model.eigenvalues
	unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))
	unstable_ωᵣ = circshift(map(x->real(x[1]),evals), (10,0))

	ωᵣlist = map(x->unstable_ωᵣ[x], 1:length(evals))
	γlist = map(x->unstable_γ[x], 1:length(evals))

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax = Axis(f[1, 1], title = "ωᵣ vs. γ", xlabel = "γ", ylabel = "ωᵣ")
	scatter!(γlist, ωᵣlist)
	hlines!(ax, [0], color=:black)
	vlines!(ax, [0], color=:black)
	
	f
end

