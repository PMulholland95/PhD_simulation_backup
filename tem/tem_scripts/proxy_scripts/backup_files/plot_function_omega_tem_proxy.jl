using DelimitedFiles, CairoMakie, LaTeXStrings

function plotOmegaTemProxy(scanlog::String, 
			   proxyoutput::String, 
			   saveplot::String,
			   n::Int64,
			   zpn::Int64;
                           resolution=(4000,2000),
                           fontsize=60,
			   linewidth=15,
			   lw2=5,
                          )

	gsc = readdlm(scanlog, skipstart=1);
	gradn = gsc[1:n,3]
	gradn = convert(Array{Float64}, gradn)
	
	g = "g"
	dat = ".dat"

	gn = []

	for i in 1:4
		push!(gn, string(proxyoutput,g,i,dat))
	end

	gn1 = vec(readdlm(gn[1]));
	gn2 = vec(readdlm(gn[2]));
	gn3 = vec(readdlm(gn[3]));
	gn4 = vec(readdlm(gn[4]));

	#=
	qpr1 = vec(readdlm("saved_data/qpr1.dat"));
	qpr2 = vec(readdlm("saved_data/qpr2.dat"));
	qpr3 = vec(readdlm("saved_data/qpr3.dat"));
	qpr4 = vec(readdlm("saved_data/qpr4.dat"));
	=#

	b = "ban"

	ban = []

	for i in 1:4
		push!(ban, string(proxyoutput,b,i,dat))
	end

	ban1 = vec(readdlm(ban[1]));
	ban2 = vec(readdlm(ban[2]));
	ban3 = vec(readdlm(ban[3]));
	ban4 = vec(readdlm(ban[4]));

	k1 = ["0.6","0.8"]
	k2 = ["0.8","1.0"]
	k3 = ["1.0","1.2"]
	k4 = ["1.2","1.4"]
	
	fig = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 = Axis(fig[1,1], title = latexstring(L"k_y \rho =",k1[zpn]), xlabel = L"a / L_n", ylabel = L"ω / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax2 = Axis(fig[1,2], title = latexstring(L"k_y \rho =",k2[zpn]), xlabel = L"a / L_n", ylabel = L"ω / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax3 = Axis(fig[2,1], title = latexstring(L"k_y \rho =",k3[zpn]), xlabel = L"a / L_n", ylabel = L"ω / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)
	ax4 = Axis(fig[2,2], title = latexstring(L"k_y \rho =",k4[zpn]), xlabel = L"a / L_n", ylabel = L"ω / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

	#lines!(ax1, gradn, qpr1, color = :green)
	scatter!(ax1, gradn, gn1, color = :blue)
	scatter!(ax1, gradn, ban1, color = :red)
	lines!(ax1, gradn, gn1, color = :blue, linewidth=linewidth)
	lines!(ax1, gradn, ban1, color = :red, linewidth=linewidth)
	hlines!(ax1, [0], color = :black, linewidth = lw2)
	vlines!(ax1, [0], color = :black, linewidth = lw2)

	#lines!(ax2, gradn, qpr2, color = :green)
	scatter!(ax2, gradn, gn2, color = :blue)
	scatter!(ax2, gradn, ban2, color = :red)
	lines!(ax2, gradn, gn2, color = :blue, linewidth=linewidth)
	lines!(ax2, gradn, ban2, color = :red, linewidth=linewidth)
	hlines!(ax2, [0], color = :black, linewidth = lw2)
	vlines!(ax2, [0], color = :black, linewidth = lw2)

	#lines!(ax3, gradn, qpr3, color = :green)
	scatter!(ax3, gradn, gn3, color = :blue)
	scatter!(ax3, gradn, ban3, color = :red)
	lines!(ax3, gradn, gn3, color = :blue, linewidth=linewidth)
	lines!(ax3, gradn, ban3, color = :red, linewidth=linewidth)
	hlines!(ax3, [0], color = :black, linewidth = lw2)
	vlines!(ax3, [0], color = :black, linewidth = lw2)

	#lines!(ax4, gradn, qpr4, color = :green)
	scatter!(ax4, gradn, gn4, color = :blue)
	scatter!(ax4, gradn, ban4, color = :red)
	lines!(ax4, gradn, gn4, color = :blue, linewidth=linewidth)
	lines!(ax4, gradn, ban4, color = :red, linewidth=linewidth)
	hlines!(ax4, [0], color = :black, linewidth = lw2)
	vlines!(ax4, [0], color = :black, linewidth = lw2)

	save(saveplot,fig)

	fig

end
