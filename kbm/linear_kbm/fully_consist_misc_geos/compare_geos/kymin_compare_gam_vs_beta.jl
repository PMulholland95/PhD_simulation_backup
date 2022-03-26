using DelimitedFiles, CairoMakie, LaTeXStrings 

sc1 = readdlm("combo_sc_dkh", skipstart=1)
sc2 = readdlm("combo_sc_dkm", skipstart=1)
sc3 = readdlm("combo_sc_dks", skipstart=1)

markersize = 80
linewidth = 20
lw2 = 10
fontsize = 120
resolution = (4000,2200)

beta1 = [0, 0.0055236, 0.0111913, 0.0170108] 
beta1 = convert(Array{Float64}, beta1);
beta1 = 100*beta1

beta2 = [0, 0.0055668, 0.0112726, 0.0170108, 0.0231552] 
beta2 = convert(Array{Float64}, beta2);
beta2 = 100*beta2

beta3 = [0, 0.0055484, 0.011235, 0.0170720, 0.0230754] 
beta3 = convert(Array{Float64}, beta3);
beta3 = 100*beta3

gam1 = sc1[:,5];
gam1 = convert(Array{Float64}, gam1);
omg1 = sc1[:,6];
omg1 = convert(Array{Float64}, omg1);  

gam2 = sc2[:,5];
gam2 = convert(Array{Float64}, gam2);
omg2 = sc2[:,6];
omg2 = convert(Array{Float64}, omg2);  

gam3 = sc3[:,5];
gam3 = convert(Array{Float64}, gam3);
omg3 = sc3[:,6];
omg3 = convert(Array{Float64}, omg3);  

n = 3

g1a = gam1[1:n:end]; 
g2a = gam1[2:n:end];
g3a = gam1[3:n:end];

g1b = gam2[1:n:end]; 
g2b = gam2[2:n:end];
g3b = gam2[3:n:end];

g1c = gam3[1:n:end]; 
g2c = gam3[2:n:end];
g3c = gam3[3:n:end];

gm1 = [g1a g2a g3a]
gm2 = [g1b g2b g3b]
gm3 = [g1c g2c g3c]

o1a = omg1[1:n:end];
o2a = omg1[2:n:end];
o3a = omg1[3:n:end];

o1b = omg2[1:n:end];
o2b = omg2[2:n:end];
o3b = omg2[3:n:end];

o1c = omg3[1:n:end];
o2c = omg3[2:n:end];
o3c = omg3[3:n:end];

om1 = [o1a o2a o3a]
om2 = [o1b o2b o3b]
om3 = [o1c o2c o3c]

xmx1=2.5
ymx1=15
ymx2=8

xax = 0:0.5:xmx1
xax = convert(Array{Float64}, xax);
yax = 0:1:ymx1
yax = convert(Array{Float64}, yax);
yax2 = -3:0.5:ymx2
yax2 = convert(Array{Float64}, yax2);

fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"β / %", ylabel = L"γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

gs1  = scatter!(beta1, gm1[:,1], color = :red, markersize = markersize, marker = '□')
gl1 = lines!(ax1, beta1, gm1[:, 1], color = :red, linewidth = linewidth)

gs2  = scatter!(beta2, gm2[:,1], color = :blue, markersize = markersize, marker = '□')
gl2 = lines!(ax1, beta2, gm2[:, 1], color = :blue, linewidth = linewidth)

gs3  = scatter!(beta3, gm3[:,1], color = :green, markersize = markersize, marker = '□')
gl3 = lines!(ax1, beta3, gm3[:, 1], color = :green, linewidth = linewidth)

hlines!(ax1, [0], color = :black, linewidth = lw2) 
vlines!(ax1, [0], color = :black, linewidth = lw2) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)


ax2 = Axis(fig[1,2], xlabel = L"β / %", ylabel = L"ωᵣ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

oms1 = scatter!(beta1, om1[:,1], color = :red, markersize = markersize, marker = '□')
oml1 = lines!(ax2, beta1, om1[:,1], color = :red, linewidth = linewidth)

oms2 = scatter!(beta2, om2[:,1], color = :blue, markersize = markersize, marker = '□')
oml2 = lines!(ax2, beta2, om2[:,1], color = :blue, linewidth = linewidth)

oms3 = scatter!(beta3, om3[:,1], color = :green, markersize = markersize, marker = '□')
oml3 = lines!(ax2, beta3, om3[:,1], color = :green, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = lw2) 
vlines!(ax2, [0], color = :black, linewidth = lw2)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

xlims!(ax1, 0, xmx1)
ylims!(ax1, -0.4, ymx1)
xlims!(ax2, 0, xmx1)
ylims!(ax2, -3, ymx2)

save("ky_compare_test.png",fig)

fig
