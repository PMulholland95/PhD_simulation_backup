using DelimitedFiles, CairoMakie, LaTeXStrings

function plotGeoCompare(proxyoutput1::String, 
			   proxyoutput2::String,
			   saveplot::String;
                           resolution=(4000,2000),
                           fontsize=60,
			   linewidth=10,
			   markersize=20,
			   lw=5)

	B1 = vec(readdlm(string(proxyoutput1,"B_z.dat")))
	B1 = convert(Array{Float64}, B1)
	κ1 = readdlm(string(proxyoutput1,"kappa_z.dat"))
	κ1 = vec(convert(Array{Float64}, κ1))
	gyy1 = vec(readdlm(string(proxyoutput1,"gyy_z.dat")))
	gyy1 = convert(Array{Float64}, gyy1)

	B2 = vec(readdlm(string(proxyoutput2,"B_z.dat")))
	B2 = convert(Array{Float64}, B2)
	κ2 = readdlm(string(proxyoutput2,"kappa_z.dat"))
	κ2 = vec(convert(Array{Float64}, κ2))
	gyy2 = vec(readdlm(string(proxyoutput2,"gyy_z.dat")))
	gyy2 = convert(Array{Float64}, gyy2)

	z = 1:1:length(B1)

##

	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(f[1,1], title = latexstring("Magnetic field strength ", L"B"), xlabel = L"z", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax2 = Axis(f[1,2], title = latexstring("Field line curvature ", L"\kappa"), xlabel = L"a / L_n", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax3 = Axis(f[2,1], title = latexstring(L"g_{yy}"), xlabel = L"z", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

	scatter!(ax1, z, B1, color = :blue, markersize=markersize)
	lines!(ax1, z, B1, color = :blue, linewidth=linewidth)
	scatter!(ax1, z, B2, color = :red, markersize=markersize)
	lines!(ax1, z, B2, color = :red, linewidth=linewidth)
	hlines!(ax1, [0], color = :black, linewidth = lw)
	vlines!(ax1, [0], color = :black, linewidth = lw)

	scatter!(ax2, z, κ1, color = :blue, markersize=markersize)
	lines!(ax2, z, κ1, color = :blue, linewidth=linewidth)
	scatter!(ax2, z, κ2, color = :red, markersize=markersize)
	lines!(ax2, z, κ2, color = :red, linewidth=linewidth)
	hlines!(ax2, [0], color = :black, linewidth = lw)
	vlines!(ax2, [0], color = :black, linewidth = lw)

	scatter!(ax3, z, gyy1, color = :blue, markersize=markersize)
	lines!(ax3, z, gyy1, color = :blue, linewidth=linewidth)
	scatter!(ax3, z, gyy2, color = :red, markersize=markersize)
	lines!(ax3, z, gyy2, color = :red, linewidth=linewidth)
	hlines!(ax3, [0], color = :black, linewidth = lw)
	vlines!(ax3, [0], color = :black, linewidth = lw)
	
	f

	save(saveplot,f)

end
