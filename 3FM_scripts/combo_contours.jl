using CairoMakie
using LaTeXStrings

function plotModeContours(model::PlasmaTurbulenceSaturationModel.FluidModel;
                           resolution=(900,600),
                           fontsize=25,
                           colormap=:vik,
                           clevels=25,
			   title= "Growth rates",
                          )
  kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
  ky = PlasmaTurbulenceSaturationModel.ky(model)
  evals = model.eigenvalues
  unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))
  unstable_ωᵣ = circshift(map(x->real(x[1]),evals), (10,0))
  marginal_ωᵣ = circshift(map(x->real(x[3]),evals), (10,0))
  
	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1, 1], title = "ωᵣ - (un)stable modes", xlabel = L"kx", ylabel = L"ky")
	cofreq = contourf!(kx, ky, unstable_ωᵣ, levels = clevels, colormap = colormap)

	Axis(f[1, 3], title = "γ - unstable modes", xlabel = L"kx", ylabel = L"ky")
	cogam = contourf!(kx, ky, unstable_γ, levels = clevels, colormap = colormap)

	Axis(f[2, 1], title = "ωᵣ - marginal modes", xlabel = L"kx", ylabel = L"ky")
	comarg = contourf!(kx, ky, marginal_ωᵣ, levels = clevels, colormap = colormap)

	Colorbar(f[1, 2], cofreq)

	Colorbar(f[1, 4], cogam)

	Colorbar(f[2, 2], comarg)
	
	f

end
