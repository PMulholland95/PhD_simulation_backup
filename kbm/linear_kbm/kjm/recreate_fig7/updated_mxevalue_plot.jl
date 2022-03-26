using DelimitedFiles, CairoMakie, LaTeXStrings 

sc = readdlm("scan.log");

markersize = 100
linewidth = 15
fontsize = 100
lw2 = 4 
resolution = (4000,1800)

n=8
c=5

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);
beta = 100*beta

beta_a = beta[1:c]
beta_b = beta[c+1:n]

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

kyrho=range(0.05,0.8, length=16); 

xax = 0:0.5:3.5
xax = convert(Array{Float64}, xax);
yax = 0:0.05:0.5
yax = convert(Array{Float64}, yax);
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

gbx = [gb1, gb2, gb3, gb4, gb5, gb6, gb7, gb8]

ob1 = omg[1:n:end];
ob2 = omg[2:n:end];
ob3 = omg[3:n:end];
ob4 = omg[4:n:end];
ob5 = omg[5:n:end];
ob6 = omg[6:n:end];
ob7 = omg[7:n:end];
ob8 = omg[8:n:end];

gmx = [maximum(gb1[.!isnan.(gb1)]),maximum(gb2[.!isnan.(gb2)]),maximum(gb3[.!isnan.(gb3)]),maximum(gb4[.!isnan.(gb4)]),maximum(gb5[.!isnan.(gb5)]),maximum(gb6[.!isnan.(gb6)]),maximum(gb7[.!isnan.(gb7)]), maximum(gb8[.!isnan.(gb8)])] 

# maxgam = convert(Array{Float64}, maxgam); 

function mxindx(x)
 	for i in 1:16
 		if gmx[x] == (gbx[x])[i]
 		return i	
 		end
 	end
end

omx = [ob1[mxindx(1)], ob2[mxindx(2)], ob3[mxindx(3)], ob4[mxindx(4)], ob5[mxindx(5)], ob6[mxindx(6)], ob7[mxindx(7)], ob8[mxindx(8)]] 

t1=20
t2=20

xmx1=3.8 
xmx2=3.8 
ymx1=0.55
ymx2=0.55

# plot(kyrho, gambet0)
fig = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"β / %", ylabel = L"γ / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = t1, yticksize = t1, xminorticksize = t2, yminorticksize = t2)

gmx_itg = scatter!(ax1, beta_a, gmx[1:c], color = :green, markersize = markersize, marker = '△') 
gmx_kbm = scatter!(ax1, beta_b, gmx[c+1:n], color = :green, markersize = markersize, marker = '□')
gmxl = lines!(ax1, beta, gmx, color = :green, linewidth = linewidth)

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth) 
hlines!(ax1, [ymx1], color = :black, linewidth = linewidth) 
vlines!(ax1, [xmx1], color = :black, linewidth = linewidth) 
hlines!(ax1, yax, color = :black, linewidth = lw2)
vlines!(ax1, xax, color = :black, linewidth = lw2)

ax2 = Axis(fig[1,2], xlabel = L"β / %", ylabel = L"ω / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = t1, yticksize = t1, xminorticksize = t2, yminorticksize = t2)

scatter!(beta_a, omx[1:c], color = :green, markersize = markersize, marker = '△')
scatter!(beta_b, omx[c+1:n], color = :green, markersize = markersize, marker = '□')
lines!(ax2, beta_a, omx[1:c], color = :green, linewidth = linewidth) 
lines!(ax2, beta_b, omx[c+1:n], color = :green, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth) 
hlines!(ax2, [ymx2], color = :black, linewidth = linewidth) 
vlines!(ax2, [xmx2], color = :black, linewidth = linewidth) 
hlines!(ax2, yax, color = :black, linewidth = lw2)
vlines!(ax2, xax, color = :black, linewidth = lw2)

#Legend(fig[1,3], [gmxal, gmxl, gmxa_itg, gmxa_kbm], [L" a/L_T = 2", L" a/L_T = 3", " ITG", " KBM"], labelsize = 40)

xlims!(ax1, 0, xmx1)
xlims!(ax2, 0, xmx2) 
ylims!(ax1, 0, ymx1) 
ylims!(ax2, 0, ymx2)

save("test2.png",fig)

fig
