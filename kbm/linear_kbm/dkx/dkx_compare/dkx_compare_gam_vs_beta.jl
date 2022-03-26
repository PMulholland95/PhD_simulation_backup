using DelimitedFiles, CairoMakie, LaTeXStrings 

sca = readdlm("sc_dkh_ky_0.025", skipstart=1)
scb = readdlm("sc_dkm_ky_0.025", skipstart=1)
scc = readdlm("sc_dks_ky_0.025", skipstart=1)

#=
scd = readdlm("sc_dkh_ky_0.8", skipstart=1)
sce = readdlm("sc_dkm_ky_0.8", skipstart=1)
scf = readdlm("sc_dks_ky_0.65", skipstart=1)
=#

markersize = 80
linewidth = 15
fontsize = 120
resolution = (4000,2200)

n=15

beta = 0.002:0.002:0.03
beta = convert(Array{Float64}, beta);
beta = 100*beta

#=
betaA = beta[1:3]
betaB = beta[4:15]
=#

gama = sca[:,7];
gama = convert(Array{Float64}, gama);
omga = sca[:,8];
omga = convert(Array{Float64}, omga);  

gamb = scb[:,7];
gamb = convert(Array{Float64}, gamb);
omgb = scb[:,8];
omgb = convert(Array{Float64}, omgb);  

gamc = scc[:,7];
gamc = convert(Array{Float64}, gamc);
omgc = scc[:,8];
omgc = convert(Array{Float64}, omgc);  

gamd = scd[:,5];
gamd = convert(Array{Float64}, gamd);
omgd = scd[:,6];
omgd = convert(Array{Float64}, omgd);  

game = sce[:,5];
game = convert(Array{Float64}, game);
omge = sce[:,6];
omge = convert(Array{Float64}, omge);  

gamf = scf[:,7];
gamf = convert(Array{Float64}, gamf);
omgf = scf[:,8];
omgf = convert(Array{Float64}, omgf);  

xmx1=3.25
ymx1=0.22
ymx2=0.05


xax = 0:0.5:xmx1
xax = convert(Array{Float64}, xax);
yax = 0:0.02:ymx1
yax = convert(Array{Float64}, yax);
yax2 = 0:0.02:ymx2
yax2 = convert(Array{Float64}, yax2);

fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"β / %", ylabel = L"γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

gma_kbm = scatter!(beta, gama, color = :red, markersize = markersize, marker = '□')
gma_kbm_l1 = lines!(ax1, beta, gama, color = :red, linewidth = linewidth)

gmb_kbm = scatter!(beta, gamb, color = :blue, markersize = markersize, marker = '□')
gmb_kbm_l1 = lines!(ax1, beta, gamb, color = :blue, linewidth = linewidth)

gmc_kbm = scatter!(beta, gamc, color = :green, markersize = markersize, marker = '□')
gmc_kbm_l1 = lines!(ax1, beta, gamc, color = :green, linewidth = linewidth)

#=
gmd_kbm = scatter!(beta, gamd, color = :red, markersize = markersize, marker = '△')
gmd_kbm_l1 = lines!(ax1, beta, gamd, color = :red, linewidth = linewidth, linestyle = :dash)

gme_kbm = scatter!(beta, game, color = :blue, markersize = markersize, marker = '△')
gme_kbm_l1 = lines!(ax1, beta, game, color = :blue, linewidth = linewidth, linestyle = :dash)

gmf_kbm = scatter!(beta, gamf, color = :green, markersize = markersize, marker = '△')
gmf_kbm_l1 = lines!(ax1, beta, gamf, color = :green, linewidth = linewidth, linestyle = :dash)
=#

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

ax2 = Axis(fig[1,2], xlabel = L"β / %", ylabel = L"ωᵣ / (c_s/a)", yticks = 0:0.01:0.045, xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

oma_kbm = scatter!(beta, omga, color = :red, markersize = markersize, marker = '□')
oma_kbm_l1 = lines!(ax2, beta, omga, color = :red, linewidth = linewidth)

omb_kbm = scatter!(beta, omgb, color = :blue, markersize = markersize, marker = '□')
omb_kbm_l1 = lines!(ax2, beta, omgb, color = :blue, linewidth = linewidth)

omc_kbm = scatter!(beta, omgc, color = :green, markersize = markersize, marker = '□')
omc_kbm_l1 = lines!(ax2, beta, omgc, color = :green, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

xlims!(ax1, 0, xmx1)
ylims!(ax1, 0, ymx1)
xlims!(ax2, 0, xmx1)
ylims!(ax2, 0, ymx2)

save("dkx_test.png",fig)

fig
