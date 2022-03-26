#using DelimitedFiles, CairoMakie, GLMakie
#GLMakie.activate!()
using CairoMakie
using LaTeXStrings

function plotGammaContours(model::PlasmaTurbulenceSaturationModel.FluidModel;
                           resolution=(3000,2000),
                           fontsize=100,
                           colormap=:vik,
                           clevels=30,
			   title= L"Growth \; rates",
                          )
  kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
  ky = PlasmaTurbulenceSaturationModel.ky(model)
  evals = model.eigenvalues
  unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))
	
  t1 = 50
  t2 = 30

  f = Figure(resolution=resolution)

  fontsize_theme = Theme(fontsize = fontsize)
  set_theme!(fontsize_theme)

  Axis(f[1, 1], title = title, xlabel = L"k_x", ylabel = L"k_y", xminorticksvisible = true, yminorticksvisible = true, xticksize = t1, yticksize = t1, xminorticksize = t2, yminorticksize = t2)
  co = contourf!(kx, ky, unstable_γ, levels = clevels, colormap = colormap)

  Colorbar(f[1, 2], co, label=L"\gamma (c_s/a)")

  save("growth_contour_3fm.png", f)

  f

end
