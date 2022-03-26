using DelimitedFiles, CairoMakie, LaTeXStrings

dsc = readdlm("d3d_scan.log");
n=11
gradn = dsc[1:n,3]
gradn = convert(Array{Float64}, gradn)

gn1 = vec(readdlm("saved_data/d1.dat"));
gn2 = vec(readdlm("saved_data/d2.dat"));
gn3 = vec(readdlm("saved_data/d3.dat"));
gn4 = vec(readdlm("saved_data/d4.dat"));

#=
qpr1 = vec(readdlm("saved_data/qpr1.dat"));
qpr2 = vec(readdlm("saved_data/qpr2.dat"));
qpr3 = vec(readdlm("saved_data/qpr3.dat"));
qpr4 = vec(readdlm("saved_data/qpr4.dat"));
=#

ban1 = vec(readdlm("saved_data/ban1.dat"));
ban2 = vec(readdlm("saved_data/ban2.dat"));
ban3 = vec(readdlm("saved_data/ban3.dat"));
ban4 = vec(readdlm("saved_data/ban4.dat"));

fig = Figure(resolution=(4000,2000))

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

#lines!(ax1, gradn, qpr1, color = :green)
scatter!(ax1, gradn, gn1, color = :blue)
scatter!(ax1, gradn, ban1, color = :red)
lines!(ax1, gradn, gn1, color = :blue)
lines!(ax1, gradn, ban1, color = :red)

#lines!(ax2, gradn, qpr2, color = :green)
scatter!(ax2, gradn, gn2, color = :blue)
scatter!(ax2, gradn, ban2, color = :red)
lines!(ax2, gradn, gn2, color = :blue)
lines!(ax2, gradn, ban2, color = :red)

#lines!(ax3, gradn, qpr3, color = :green)
scatter!(ax3, gradn, gn3, color = :blue)
scatter!(ax3, gradn, ban3, color = :red)
lines!(ax3, gradn, gn3, color = :blue)
lines!(ax3, gradn, ban3, color = :red)

#lines!(ax4, gradn, qpr4, color = :green)
scatter!(ax4, gradn, gn4, color = :blue)
scatter!(ax4, gradn, ban4, color = :red)
lines!(ax4, gradn, gn4, color = :blue)
lines!(ax4, gradn, ban4, color = :red)

save("phi_proxy_test.png",fig)

fig
