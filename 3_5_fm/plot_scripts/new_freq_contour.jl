using CairoMakie
using LaTeXStrings

function plotFreqContours(model::PlasmaTurbulenceSaturationModel.FluidModel;
                           resolution=(3000,2000),
                           fontsize=100,
                           colormap=:vik,
                           clevels=30,
			   title= L"(Un)stable \; mode \; frequency",
                          )
  kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
  ky = PlasmaTurbulenceSaturationModel.ky(model)
  evals = model.eigenvalues
  unstable_ωᵣ = circshift(map(x->real(x[1]),evals), (10,0))
  
	t1 = 50
	t2 = 30

	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1, 1], 
	     title = title, 
	     xlabel = L"kx", 
	     ylabel = L"ky", 
	     xminorticksvisible = true, 
	     yminorticksvisible = true, 
	     xticksize = t1, 
	     yticksize = t1, 
	     xminorticksize = t2, 
	     yminorticksize = t2)

	cofreq = contourf!(kx, ky, unstable_ωᵣ, levels = clevels, colormap = colormap)

	Colorbar(f[1, 2], cofreq, label=L"\omega (c_s/a)")

	save("freq_contour_3fm.png",f)
	
	f

end
