using DelimitedFiles, CairoMakie, LaTeXStrings 

sc = readdlm("combo_sc", skipstart=1)

markersize = 80
linewidth = 20
lw2 = 10
fontsize = 120
resolution = (4000,2200)

beta = [0, 0.0056383, 0.0114151, 0.0173420, 0.0234322] 
beta = convert(Array{Float64}, beta);
beta = 100*beta

gam = sc[:,5];
gam = convert(Array{Float64}, gam);
omg = sc[:,6];
omg = convert(Array{Float64}, omg);  

n = 3

g1 = gam[1:n:end]; 
g2 = gam[2:n:end];
g3 = gam[3:n:end];

gm = [g1 g2 g3]

o1 = omg[1:n:end];
o2 = omg[2:n:end];
o3 = omg[3:n:end];

om = [o1 o2 o3]

xmx1=2.5
ymx1=7.5
ymx2=7

xax = 0:0.5:xmx1
xax = convert(Array{Float64}, xax);
yax = 0:0.5:ymx1
yax = convert(Array{Float64}, yax);
yax2 = 0:0.5:ymx2
yax2 = convert(Array{Float64}, yax2);

fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"β / %", ylabel = L"γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

gs1  = scatter!(beta, gm[:,1], color = :red, markersize = markersize, marker = '□')
gl1 = lines!(ax1, beta, gm[:, 1], color = :red, linewidth = linewidth)

gs2  = scatter!(beta, gm[:,2], color = :blue, markersize = markersize, marker = '□')
gl2 = lines!(ax1, beta, gm[:, 2], color = :blue, linewidth = linewidth)

gs3  = scatter!(beta, gm[:,3], color = :green, markersize = markersize, marker = '□')
gl3 = lines!(ax1, beta, gm[:, 3], color = :green, linewidth = linewidth)

hlines!(ax1, [0], color = :black, linewidth = lw2) 
vlines!(ax1, [0], color = :black, linewidth = lw2) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

ax2 = Axis(fig[1,2], xlabel = L"β / %", ylabel = L"ωᵣ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

oms1 = scatter!(beta, om[:,1], color = :red, markersize = markersize, marker = '□')
oml1 = lines!(ax2, beta, om[:,1], color = :red, linewidth = linewidth)

oms2 = scatter!(beta, om[:,2], color = :blue, markersize = markersize, marker = '□')
oml2 = lines!(ax2, beta, om[:,2], color = :blue, linewidth = linewidth)

oms3 = scatter!(beta, om[:,3], color = :green, markersize = markersize, marker = '□')
oml3 = lines!(ax2, beta, om[:,3], color = :green, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = lw2) 
vlines!(ax2, [0], color = :black, linewidth = lw2)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

xlims!(ax1, 0, xmx1)
ylims!(ax1, 0, ymx1)
xlims!(ax2, 0, xmx1)
ylims!(ax2, 0, ymx2)

save("kjm_ky_compare_test.png",fig)

fig
