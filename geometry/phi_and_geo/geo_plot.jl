using DelimitedFiles, CairoMakie, LaTeXStrings

lw = 15
fs = 120
resolution=(3000,2000)

sc = readdlm("../../new_norm_geos/ksenia_norm_gist_incl_beta/fis_gist/gist_tube_wout_FIS_4_norm.txt_ksu", skipstart=15)

B = sc[:,4]
B = convert(Array{Float64, 1}, B);
#maxB = maximum(B)
#B = B/maxB

κ = sc[:,6]
κ = convert(Array{Float64, 1}, κ);
#maxκ = maximum(κ)
#κ = κ/maxκ

nz=range(0,120, length=120)

fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fs)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], yticklabelcolor = :blue, ylabel = L"B", xlabel = L"z")
ax2 = Axis(fig[1,1], yticklabelcolor = :red, yaxisposition = :right, ylabel = L"\kappa")

lines!(ax1, nz, B, color = :blue, linewidth = lw)
lines!(ax2, nz, κ, color = :red, linewidth = lw)
hlines!(ax2, [0], color = :black, linestyle = :dash, linewidth = 10)
#lines!(ax3, nz, ϕ, color = :purple)

save("misc_data_plots/ksenia_norm_geo_plots/geo_norm_s0_0.36_fis04_ksenia.png",fig)

fig
