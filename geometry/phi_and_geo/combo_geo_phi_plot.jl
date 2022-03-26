using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

sc = readdlm("gist20");
sc1 = readdlm("edit_phi_z");

B = sc[:,4]
B = convert(Array{Float64}, B);
a = 1:1:128
C = zeros(0)

for i in a
	append!(C, B[i]-1.0)
end
B = C
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

hlines!(ax1, [0], xmax=[1], color=:black)
lines!(ax1, nz, B, color = :blue)
lines!(ax1, nz, κ, color = :red)
lines!(ax1, nz, ϕ/5, color = :purple)

ax1.title = "B, κ and ϕ [kyρ = 0.05, β = 3.5%]"

ax1.xlabel = "z"

xlims!(ax1, 0, 128)

ylims!(ax1, -0.1, 0.3)

save("combo_phi_geo_plots.png",fig)

fig
