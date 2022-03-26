using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

fontsize = 100
axsize = 100
markersize = 70
linewidth = 8

sc0 = readdlm("../scanlogs_broad/p21sc2_scan")
sca = readdlm("../scanlogs_broad/p24_sc0_scan")
scb = readdlm("../scanlogs_broad/p24_sc1_scan")

# Only use the following for p24_sc0_scan and p24_sc1_scan
#=
gama = sca[:,2];
gama = convert(Array{Float64}, gama);
omga = sca[:,3];
omga = convert(Array{Float64}, omga);  

gamb = scb[:,7];
gamb = convert(Array{Float64}, gamb);
omgb = scb[:,8];
omgb = convert(Array{Float64}, omgb);  

gam0 = [gama; gamb]
omg0 = [omga; omgb]
=#
# end of alternative approach for p24

n=15
c=13

beta = 0.002:0.002:0.03
beta = convert(Array{Float64}, beta);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);

ky = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

ky1 = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14]
ky2 = [0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

xax = 0:0.1:0.9     
xax = convert(Array{Float64}, xax);
yax = 0:0.1:0.8
yax = convert(Array{Float64}, yax);
yax1 = 0:0.01:0.12
yax1 = convert(Array{Float64}, yax1);

gb10a1 = (gam0[10:n:end])[1:c];
gb10a2 = (gam0[10:n:end])[c+1:27];
ob10a1 = (omg0[10:n:end])[1:c];
ob10a2 = (omg0[10:n:end])[c+1:27];

resolution=(4000,2000)

f = Figure(resolution=resolution)
fontsize_theme = Theme(fontsize = fontsize) 
set_theme!(fontsize_theme)

ax1 = Axis(f[1,1], xlabel = L"k_y\rho", xlabelsize = axsize, ylabel = L"ω / (c_s/a)", ylabelsize = axsize, xminorticksvisible = true, yminorticksvisible = true)
scatter!(ky1, ob10a1, color = :blue, markersize = markersize, marker = '□') 
lines!(ky1, ob10a1, color = :blue, linewidth = linewidth)
scatter!(ky2, ob10a2, color = :blue, markersize = markersize, marker = '△') 
lines!(ky2, ob10a2, color = :blue, linewidth = linewidth)
hlines!(ax1, [0], color = :black, linewidth = 14)
vlines!(ax1, [0], color = :black, linewidth = 14)
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

ax2 = Axis(f[1,2], xlabel = L"k_y\rho", xlabelsize = axsize, ylabel = L"γ / (c_s/a)", ylabelsize = axsize, xminorticksvisible = true, yminorticksvisible = true)
scatter!(ky1, gb10a1, color = :red, markersize = markersize, marker = '□')
lines!(ky1, gb10a1, color = :red, linewidth = linewidth) 
scatter!(ky2, gb10a2, color = :red, markersize = markersize, marker = '△')
lines!(ky2, gb10a2, color = :red, linewidth = linewidth) 
hlines!(ax2, [0], color = :black, linewidth = 14)
vlines!(ax2, [0], color = :black, linewidth = 14)
hlines!(ax2, yax1, color = :grey)  
vlines!(ax2, xax, color = :grey) 

xlims!(ax1, 0, 0.95)
xlims!(ax2, 0, 0.95)
ylims!(ax1, 0, 0.82) 
ylims!(ax2, 0, 0.12) 

save("lin_vs_ky.png", f)

f
