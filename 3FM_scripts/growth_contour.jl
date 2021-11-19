#using DelimitedFiles, CairoMakie, GLMakie
#GLMakie.activate!()
using CairoMakie
using LaTeXStrings

function plotGammaContours(model::PlasmaTurbulenceSaturationModel.FluidModel;
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


  f = Figure(resolution=resolution)

  fontsize_theme = Theme(fontsize = fontsize)
  set_theme!(fontsize_theme)

  Axis(f[1, 1], title = title, xlabel = L"k_x", ylabel = L"k_y")
  co = contourf!(kx, ky, unstable_γ, levels = clevels, colormap = colormap)

  Colorbar(f[1, 2], co, label=L"\gamma (c_s/a)")

  f

end
