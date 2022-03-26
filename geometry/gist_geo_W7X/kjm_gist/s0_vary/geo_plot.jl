using DelimitedFiles, CairoMakie, GLMakie, LaTeXStrings
GLMakie.activate!()

lw = 15
fs = 120

sc = readdlm("kjm_0_s0_0.5.dat")

B = sc[:,4]
B = convert(Array{Float64, 1}, B);
#maxB = maximum(B)
#B = B/maxB

κ = sc[:,6]
κ = convert(Array{Float64, 1}, κ);
#maxκ = maximum(κ)
#κ = κ/maxκ

nz=range(0,128, length=128)

fig = Figure()

fontsize_theme = Theme(fontsize = fs)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], yticklabelcolor = :blue, ylabel = L"B", ylabelcolor = :blue, xlabel = L"z")
ax2 = Axis(fig[1,1], yticklabelcolor = :red, yaxisposition = :right, ylabel = L"\kappa", ylabelcolor = :red)

lines!(ax1, nz, B, color = :blue, linewidth = lw)
lines!(ax2, nz, κ, color = :red, linewidth = lw)
hlines!(ax2, [0], color = :black, linestyle = :dash, linewidth = 10)
#lines!(ax3, nz, ϕ, color = :purple)

save("geo_plots.png",fig)

fig
