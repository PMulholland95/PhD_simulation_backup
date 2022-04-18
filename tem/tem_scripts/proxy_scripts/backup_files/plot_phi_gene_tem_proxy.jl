using DelimitedFiles, CairoMakie, LaTeXStrings

linewidth = 15

gsc = readdlm("start_data/d3d_scan.log");
n=11
gradn = dsc[1:n,3]
gradn = convert(Array{Float64}, gradn)

gn1 = vec(readdlm("saved_data/save_phi_gene/gd1.dat"));
gn2 = vec(readdlm("saved_data/save_phi_gene/gd2.dat"));
gn3 = vec(readdlm("saved_data/save_phi_gene/gd3.dat"));
gn4 = vec(readdlm("saved_data/save_phi_gene/gd4.dat"));

#=
qpr1 = vec(readdlm("saved_data/qpr1.dat"));
qpr2 = vec(readdlm("saved_data/qpr2.dat"));
qpr3 = vec(readdlm("saved_data/qpr3.dat"));
qpr4 = vec(readdlm("saved_data/qpr4.dat"));
=#

ban1 = vec(readdlm("saved_data/save_phi_gene/gban1.dat"));
ban2 = vec(readdlm("saved_data/save_phi_gene/gban2.dat"));
ban3 = vec(readdlm("saved_data/save_phi_gene/gban3.dat"));
ban4 = vec(readdlm("saved_data/save_phi_gene/gban4.dat"));

fig = Figure(resolution=(4000,2000))

fontsize_theme = Theme(fontsize = 50)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

#lines!(ax1, gradn, qpr1, color = :green)
scatter!(ax1, gradn, gn1, color = :blue, linewidth=linewidth)
scatter!(ax1, gradn, ban1, color = :red, linewidth=linewidth)
lines!(ax1, gradn, gn1, color = :blue, linewidth=linewidth)
lines!(ax1, gradn, ban1, color = :red, linewidth=linewidth)

#lines!(ax2, gradn, qpr2, color = :green)
scatter!(ax2, gradn, gn2, color = :blue, linewidth=linewidth)
scatter!(ax2, gradn, ban2, color = :red, linewidth=linewidth)
lines!(ax2, gradn, gn2, color = :blue, linewidth=linewidth)
lines!(ax2, gradn, ban2, color = :red, linewidth=linewidth)

#lines!(ax3, gradn, qpr3, color = :green)
scatter!(ax3, gradn, gn3, color = :blue, linewidth=linewidth)
scatter!(ax3, gradn, ban3, color = :red, linewidth=linewidth)
lines!(ax3, gradn, gn3, color = :blue, linewidth=linewidth)
lines!(ax3, gradn, ban3, color = :red, linewidth=linewidth)

#lines!(ax4, gradn, qpr4, color = :green)
scatter!(ax4, gradn, gn4, color = :blue, linewidth=linewidth)
scatter!(ax4, gradn, ban4, color = :red, linewidth=linewidth)
lines!(ax4, gradn, gn4, color = :blue, linewidth=linewidth)
lines!(ax4, gradn, ban4, color = :red, linewidth=linewidth)

save("phi_gene_proxy_test.png",fig)

fig
