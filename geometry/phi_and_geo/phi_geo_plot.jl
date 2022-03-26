using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

sc = readdlm("gist20");
sc1 = readdlm("edit_phi_z");

B = sc[:,4]
B = convert(Array{Float64}, B);
#maxB = maximum(B)
#B = B/maxB

κ = sc[:,6]
κ = convert(Array{Float64}, κ);
#maxκ = maximum(κ)
#κ = κ/maxκ

ϕ = sc1[:,3]
ϕ = convert(Array{Float64}, ϕ);
maxϕ = maximum(ϕ)
ϕ = ϕ/maxϕ


nz=range(0,128, length=128)

fig = Figure()

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[1,3])

lines!(ax1, nz, B, color = :blue)

hlines!(ax2, [0], xmax=[1], color=:black)
lines!(ax2, nz, κ, color = :red)

lines!(ax3, nz, ϕ, color = :purple)

ax1.title = "B"
ax2.title = "κ"
ax3.title = "ϕ [kyρ = 0.05, β = 3.5%]"

ax1.xlabel = "z"
ax2.xlabel = "z"
ax3.xlabel = "z"

xlims!(ax1, 0, 128)
xlims!(ax2, 0, 128)
xlims!(ax2, 0, 128)

ylims!(ax1, 0.85, 1.3)
ylims!(ax2, -0.1, 0.3)
ylims!(ax3, 0.0, 1.1)

save("phi_geo_plots.png",fig)

fig
