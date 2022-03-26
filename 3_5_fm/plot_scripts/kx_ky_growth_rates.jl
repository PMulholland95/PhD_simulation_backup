using CairoMakie
using LaTeXStrings

function plotGammaWavenumber(model::PlasmaTurbulenceSaturationModel.FluidModel;
		kx=Float64,
		ky=Float64,
		resolution=(5000,8000),
		res2=(3000,2000),
		fontsize=100,
		title1=L"Growth \; rate \; \gamma (kx)",
		title2=L"Growth \; rate \; \gamma (ky)",
		)
	
	nkx, nky = PTSM.index(kx,ky,model.spectralGrid)
	kx = PTSM.kx(model,negativeModes=true)
 	ky = PTSM.ky(model)

	evals = model.eigenvalues
	unstable_γ = map(x->imag(x[1]),evals)
	xshift_unstable_γ = circshift(map(x->imag(x[1]),evals), (10,0))

	gkx = xshift_unstable_γ[1:length(kx),nky] 
	gky = unstable_γ[nkx,1:length(ky)] 

	lw=15
	lw2=2
	xax = -0.4:0.1:0.4
	xax = convert(Array{Float64}, xax);
	yax = 0.20:0.005:0.215
	yax = convert(Array{Float64}, yax);
	xax2 = -0:0.1:0.5
	xax2 = convert(Array{Float64}, xax2);
	yax2 = 0:0.05:0.25
	yax2 = convert(Array{Float64}, yax2);
	xmx1 = 0.52
	ymx1 = 0.216
	xmx2 = 0.52
	ymx2 = 0.22
	t1 = 50
	t2 = 30
	
	# Combined figures
	
	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(f[1,1], title = title1, xlabel = L"k_x")
	lines!(ax1, kx, gkx, color = :red, linewidth = lw)
	hlines!(ax1, [0.198], color = :black, linewidth = lw)
	vlines!(ax1, [0], color = :black, linewidth = lw2)
	hlines!(ax1, [ymx1], color = :black, linewidth = lw)
	vlines!(ax1, [xmx1], color = :black, linewidth = lw)
	vlines!(ax1, [-xmx1], color = :black, linewidth = lw)
	hlines!(ax1, yax, color = :black, linewidth = lw2)
	vlines!(ax1, xax, color = :black, linewidth = lw2)

	ax2 = Axis(f[2,1], title = title2, xlabel = L"k_y")  
	lines!(ax2, ky, gky, color = :blue, linewidth = lw)
	hlines!(ax2, [0], color = :black, linewidth = lw)
	vlines!(ax2, [0], color = :black, linewidth = lw)
	hlines!(ax2, [ymx2], color = :black, linewidth = lw)
	vlines!(ax2, [xmx2], color = :black, linewidth = lw)
	hlines!(ax2, yax2, color = :black, linewidth = lw2)
	vlines!(ax2, xax2, color = :black, linewidth = lw2)

	
	xlims!(ax1, -xmx1, xmx1)
	xlims!(ax2, 0, xmx2)
	ylims!(ax1, 0.198, ymx1)
	ylims!(ax2, 0, ymx2)


	save("kx_ky_growth_rates.png",f)

	f

	# Now inividual figures

	f1 = Figure(resolution=res2)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(f1[1,1], 
		   xlabel = L"k_x", 
		   ylabel = L"\gamma / (c_s/a)", 
		   xminorticksvisible = true, 
		   yminorticksvisible = true, 
		   xticksize = t1, 
		   yticksize = t1, 
		   xminorticksize = t2, 
		   yminorticksize = t2)

	lines!(ax1, kx, gkx, color = :red, linewidth = lw)
	hlines!(ax1, [0.198], color = :black, linewidth = lw)
	vlines!(ax1, [0], color = :black, linewidth = lw2)
	hlines!(ax1, [ymx1], color = :black, linewidth = lw)
	vlines!(ax1, [xmx1], color = :black, linewidth = lw)
	vlines!(ax1, [-xmx1], color = :black, linewidth = lw)
	hlines!(ax1, yax, color = :black, linewidth = lw2)
	vlines!(ax1, xax, color = :black, linewidth = lw2)

	xlims!(ax1, -xmx1, xmx1)
	ylims!(ax1, 0.198, ymx1)


	save("kx_growth_rate.png",f1)

	f1


	f2 = Figure(resolution=res2)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax2 = Axis(f2[1,1], 
		   xlabel = L"k_y", 
		   ylabel = L"\gamma / (c_s/a)", 
		   xminorticksvisible = true, 
		   yminorticksvisible = true, 
		   xticksize = t1, 
		   yticksize = t1, 
		   xminorticksize = t2, 
		   yminorticksize = t2)
  
	lines!(ax2, ky, gky, color = :blue, linewidth = lw)
	hlines!(ax2, [0], color = :black, linewidth = lw)
	vlines!(ax2, [0], color = :black, linewidth = lw)
	hlines!(ax2, [ymx2], color = :black, linewidth = lw)
	vlines!(ax2, [xmx2], color = :black, linewidth = lw)
	hlines!(ax2, yax2, color = :black, linewidth = lw2)
	vlines!(ax2, xax2, color = :black, linewidth = lw2)

	
	xlims!(ax2, 0, xmx2)
	ylims!(ax2, 0, ymx2)


	save("ky_growth_rate.png",f2)

	f2

end


