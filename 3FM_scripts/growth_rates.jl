using CairoMakie
using LaTeXStrings

function plotGammaWavenumber(model::PlasmaTurbulenceSaturationModel.FluidModel;
		nkx=Int,
		nky=Int,
		resolution=(1600,900),
		fontsize=25,
		title="Growth rates",
		)
	
	kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
 	ky = PlasmaTurbulenceSaturationModel.ky(model)
	evals = model.eigenvalues
	unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))

	gkx = unstable_γ[1:length(kx),nky] 
	gky = unstable_γ[nkx,1:length(ky)] 

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1,1], title = title, xlabel = L"k_x")
	lines!(kx, gkx, color = :red)

	Axis(f[1,2], title = title, xlabel = L"k_y")  
	lines!(ky, gky, color = :blue)

	f
end

