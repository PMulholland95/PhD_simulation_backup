using CairoMakie
using LaTeXStrings

function plotMargContours(model::PlasmaTurbulenceSaturationModel.FluidModel;
                           resolution=(3000,2000),
                           fontsize=100,
                           colormap=:vik,
                           clevels=30,
			   title= L"Marginal \; mode \; frequency",
                          )
  kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
  ky = PlasmaTurbulenceSaturationModel.ky(model)
  evals = model.eigenvalues
  marginal_ωᵣ = circshift(map(x->real(x[3]),evals), (10,0))
  
	t1 = 50
	t2 = 30

	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1, 1], title = title, xlabel = L"kx", ylabel = L"ky", xminorticksvisible = true, yminorticksvisible = true, xticksize = t1, yticksize = t1, xminorticksize = t2, yminorticksize = t2)
	comarg = contourf!(kx, ky, marginal_ωᵣ, levels = clevels, colormap = colormap)

	Colorbar(f[1, 2], comarg, label=L"\omega (c_s/a)")

	save("marg_contour_3fm.png",f)
	
	f

end
