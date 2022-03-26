using DelimitedFiles, CairoMakie, GLMakie, LaTeXStrings 
GLMakie.activate!()

sc = readdlm("../scanlogs_broad/p01sc0_scan");

markersize = 80
linewidth = 12
fontsize = 80
resolution= (3000,2000)

n=15

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);
beta = 100*beta

beta6a = beta[1:6]
beta6b = beta[7:15]

gam = sc[:,7];
gam = convert(Array{Float64}, gam);
omg = sc[:,8];
omg = convert(Array{Float64}, omg);

kyrho= 0.05:0.05:0.9;

xax = 0:0.5:3
xax = convert(Array{Float64}, xax);
yax = 0:0.05:0.4
yax = convert(Array{Float64}, yax);
yax2 = 0:0.05:0.6
yax2 = convert(Array{Float64}, yax2);
# We now want to gather together gamma values for range of kyrho, and fixed beta
# In the end, only want to plot: gamma|max vs. beta 

gb1 = gam[1:n:end];
gb2 = gam[2:n:end];
gb3 = gam[3:n:end];
gb4 = gam[4:n:end];
gb5 = gam[5:n:end];
gb6 = gam[6:n:end];
gb7 = gam[7:n:end];
gb8 = gam[8:n:end];
gb9 = gam[9:n:end];
gb10 = gam[10:n:end];
gb11 = gam[11:n:end];
gb12 = gam[12:n:end];
gb13 = gam[13:n:end];
gb14 = gam[14:n:end];
gb15 = gam[15:n:end];

gbx = [gb1, gb2, gb3, gb4, gb5, gb6, gb7, gb8, gb9, gb10, gb11, gb12, gb13, gb14, gb15]

ob1 = omg[1:n:end];
ob2 = omg[2:n:end];
ob3 = omg[3:n:end];
ob4 = omg[4:n:end];
ob5 = omg[5:n:end];
ob6 = omg[6:n:end];
ob7 = omg[7:n:end];
ob8 = omg[8:n:end];
ob9 = omg[9:n:end];
ob10 = omg[10:n:end];
ob11 = omg[11:n:end];
ob12 = omg[12:n:end];
ob13 = omg[13:n:end];
ob14 = omg[14:n:end];
ob15 = omg[15:n:end];

gmx = [maximum(gb1[.!isnan.(gb1)]),maximum(gb2[.!isnan.(gb2)]),maximum(gb3[.!isnan.(gb3)]),maximum(gb4[.!isnan.(gb4)]),maximum(gb5[.!isnan.(gb5)]),maximum(gb6[.!isnan.(gb6)]),maximum(gb7[.!isnan.(gb7)]),maximum(gb8[.!isnan.(gb8)]),maximum(gb9[.!isnan.(gb9)]),maximum(gb10[.!isnan.(gb10)]),maximum(gb11[.!isnan.(gb11)]),maximum(gb12[.!isnan.(gb12)]),maximum(gb13[.!isnan.(gb13)]),maximum(gb14[.!isnan.(gb14)]),maximum(gb15[.!isnan.(gb15)]),] 
# maxgam = convert(Array{Float64}, maxgam); 

function mxindx(x)
 	for i in 1:28
 		if gmx[x] == (gbx[x])[i]
 		return i	
 		end
 	end
end

omx = [ob1[mxindx(1)], ob2[mxindx(2)], ob3[mxindx(3)], ob4[mxindx(4)], ob5[mxindx(5)], ob6[mxindx(6)], ob7[mxindx(7)], ob8[mxindx(8)], ob9[mxindx(9)], ob10[mxindx(10)], ob11[mxindx(11)], ob12[mxindx(12)], ob13[mxindx(13)], ob14[mxindx(14)], ob15[mxindx(15)]]

# plot(kyrho, gambet0)
fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"β / %", ylabel = L"γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

gmx_itg = scatter!(beta6a, gmx[1:6], color = :red, markersize = markersize, marker = '△') 
gmx_kbm = scatter!(beta6b, gmx[7:15], color = :red, markersize = markersize, marker = '□') 
gmxl1 = lines!(ax1, beta6a, gmx[1:6], color = :red, linewidth = linewidth)
gmxl2 = lines!(ax1, beta6b, gmx[7:15], color = :red, linewidth = linewidth)

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

ax2 = Axis(fig[1,2], xlabel = L"β / %", ylabel = L"ωᵣ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

scatter!(beta6a, omx[1:6], color = :red, markersize = markersize, marker = '△')
scatter!(beta6b, omx[7:15], color = :red, markersize = markersize, marker = '□')
lines!(ax2, beta6a, omx[1:6], color = :red, linewidth = linewidth) 
lines!(ax2, beta6b, omx[7:15], color = :red, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

xlims!(ax1, 0, 3.2)
xlims!(ax2, 0, 3.2)
ylims!(ax1, 0, 0.22)
ylims!(ax2, 0, 0.45)

save("gamVSbeta.png", fig)
#Legend(fig[1,3], [gmxal, gmxl, gmxa_itg, gmxa_kbm], [L" a/L_T = 2", L" a/L_T = 3", " ITG", " KBM"], labelsize = 40)

fig
