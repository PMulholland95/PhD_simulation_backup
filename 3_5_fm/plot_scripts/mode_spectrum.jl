using CairoMakie
using LaTeXStrings

xax = 0:0.05:0.25
xax = convert(Array{Float64}, xax);
yax = -0.4:0.1:0
yax = convert(Array{Float64}, yax);
xmx1 = 0.26
ymx1 = -0.42
t1 = 50
t2 = 30
mkrs = 30
lw = 12
lw2 = 3

function plotModeSpectrum(model::PlasmaTurbulenceSaturationModel.FluidModel;
		resolution=(3000,2000),
		fontsize=100,
		title="Mode spectrum",
		)

	kx = PlasmaTurbulenceSaturationModel.kx(model,negativeModes=true)
 	ky = PlasmaTurbulenceSaturationModel.ky(model)

	evals = model.eigenvalues
	unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))
	unstable_ωᵣ = circshift(map(x->real(x[1]),evals), (10,0))

	ωᵣlist = map(x->unstable_ωᵣ[x], 1:length(evals))
	γlist = map(x->unstable_γ[x], 1:length(evals))

	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(f[1, 1], 
		   xlabel = L"\gamma (c_s/a)", 
		   ylabel = L"\omega (c_s/a)",
		xminorticksvisible = true, 
		yminorticksvisible = true, 
		xticksize = t1, 
		yticksize = t1, 
		xminorticksize = t2, 
		yminorticksize = t2)

	scatter!(γlist, ωᵣlist, color=:blue, markersize = mkrs)
	hlines!(ax1, [0], color=:black, linewidth = lw)
	vlines!(ax1, [0], color=:black, linewidth = lw)
	hlines!(ax1, [0.02], color = :black, linewidth = lw)
	hlines!(ax1, [ymx1], color = :black, linewidth = lw)
	vlines!(ax1, [xmx1], color = :black, linewidth = lw)
	vlines!(ax1, [-0.02], color = :black, linewidth = lw)
	hlines!(ax1, yax, color = :black, linewidth = lw2)
	vlines!(ax1, xax, color = :black, linewidth = lw2)
	
	xlims!(ax1, -0.02, xmx1)
	ylims!(ax1, ymx1, 0.02)

	save("mode_scatter_spectrum.png", f)

	f
end

