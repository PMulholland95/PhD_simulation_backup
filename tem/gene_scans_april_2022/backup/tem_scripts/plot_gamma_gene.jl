using DelimitedFiles, CairoMakie, LaTeXStrings

function plotGamma(scanlog::String, 
		   saveplot::String; 
		   n = 11,
                   resolution=(4000,2000),
                   fontsize=60,
		   linewidth=15,
		   lw2=5,
                   )

	gsc = readdlm(scanlog, skipstart=1);
	gradn = gsc[1:n,3]
	gradn = convert(Array{Float64}, gradn)
	
	gam1 = gsc[1:n,7];
	gam2 = gsc[n+1:2n,7];
	gam3 = gsc[2n+1:3n,7];
	gam4 = gsc[3n+1:4n,7];

	gam1 = convert(Array{Float64}, gam1)
	gam2 = convert(Array{Float64}, gam2)
	gam3 = convert(Array{Float64}, gam3)
	gam4 = convert(Array{Float64}, gam4)

	omg1 = gsc[1:n,8];
	omg2 = gsc[n+1:2n,8];
	omg3 = gsc[2n+1:3n,8];
	omg4 = gsc[3n+1:4n,8];

	omg1 = convert(Array{Float64}, omg1)
	omg2 = convert(Array{Float64}, omg2)
	omg3 = convert(Array{Float64}, omg3)
	omg4 = convert(Array{Float64}, omg4)

	fig = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(fig[1,1], title = L"k_y \rho = 0.8", xlabel = L"a / L_n", ylabel = L"ω,γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax2 = Axis(fig[1,2], title = L"k_y \rho = 1.0", xlabel = L"a / L_n", ylabel = L"ω,γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax3 = Axis(fig[2,1], title = L"k_y \rho = 1.2", xlabel = L"a / L_n", ylabel = L"ω,γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax4 = Axis(fig[2,2], title = L"k_y \rho = 1.4", xlabel = L"a / L_n", ylabel = L"ω,γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

	#lines!(ax1, gradn, qpr1, color = :green)
	scatter!(ax1, gradn, gam1, color = :green)
	scatter!(ax1, gradn, omg1, color = :blue)
	lines!(ax1, gradn, gam1, color = :green, linewidth=linewidth)
	lines!(ax1, gradn, omg1, color = :blue, linewidth=linewidth)
	hlines!(ax1, [0], color = :black, linewidth = lw2)
	vlines!(ax1, [0], color = :black, linewidth = lw2)

	#lines!(ax2, gradn, qpr2, color = :green)
	scatter!(ax2, gradn, gam2, color = :green)
	scatter!(ax2, gradn, omg2, color = :blue)
	lines!(ax2, gradn, gam2, color = :green, linewidth=linewidth)
	lines!(ax2, gradn, omg2, color = :blue, linewidth=linewidth)
	hlines!(ax2, [0], color = :black, linewidth = lw2)
	vlines!(ax2, [0], color = :black, linewidth = lw2)

	#lines!(ax3, gradn, qpr3, color = :green)
	scatter!(ax3, gradn, gam3, color = :green)
	scatter!(ax3, gradn, omg3, color = :blue)
	lines!(ax3, gradn, gam3, color = :green, linewidth=linewidth)
	lines!(ax3, gradn, omg3, color = :blue, linewidth=linewidth)
	hlines!(ax3, [0], color = :black, linewidth = lw2)
	vlines!(ax3, [0], color = :black, linewidth = lw2)

	#lines!(ax4, gradn, qpr4, color = :green)
	scatter!(ax4, gradn, gam4, color = :green)
	scatter!(ax4, gradn, omg4, color = :blue)
	lines!(ax4, gradn, gam4, color = :green, linewidth=linewidth)
	lines!(ax4, gradn, omg4, color = :blue, linewidth=linewidth)
	hlines!(ax4, [0], color = :black, linewidth = lw2)
	vlines!(ax4, [0], color = :black, linewidth = lw2)

	save(saveplot,fig)

	fig

end
