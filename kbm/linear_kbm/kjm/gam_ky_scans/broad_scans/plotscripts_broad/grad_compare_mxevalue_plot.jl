using DelimitedFiles, CairoMakie, GLMakie, LaTeXStrings 
GLMakie.activate!()


sc = readdlm("p20sc6_scan");
sc1 = readdlm("p21sc2_scan");
#include this in rest of script
sc2 = readdlm("p20sc7_scan")

markersize = 80
linewidth = 15
fontsize = 120

n=15

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);
beta = 100*beta

beta6a = beta[1:6]
beta6b = beta[7:15]
beta7a = beta[1:7]
beta7b = beta[8:15]
beta8a = beta[1:8]
beta8b = beta[9:15]

gam = sc[:,7];
gam = convert(Array{Float64}, gam);
omg = sc[:,8];
omg = convert(Array{Float64}, omg);

gam1 = sc1[:,7];
gam1 = convert(Array{Float64}, gam1);
omg1 = sc1[:,8];
omg1 = convert(Array{Float64}, omg1);

gam2 = sc2[:,7];
gam2 = convert(Array{Float64}, gam2);
omg2 = sc2[:,8];
omg2 = convert(Array{Float64}, omg2);

kyrho= [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2,0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];


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

gb1a = gam1[1:n:end];
gb2a = gam1[2:n:end];
gb3a = gam1[3:n:end];
gb4a = gam1[4:n:end];
gb5a = gam1[5:n:end];
gb6a = gam1[6:n:end];
gb7a = gam1[7:n:end];
gb8a = gam1[8:n:end];
gb9a = gam1[9:n:end];
gb10a = gam1[10:n:end];
gb11a = gam1[11:n:end];
gb12a = gam1[12:n:end];
gb13a = gam1[13:n:end];
gb14a = gam1[14:n:end];
gb15a = gam1[15:n:end];

gbxa = [gb1a, gb2a, gb3a, gb4a, gb5a, gb6a, gb7a, gb8a, gb9a, gb10a, gb11a, gb12a, gb13a, gb14a, gb15a]

ob1a = omg1[1:n:end];
ob2a = omg1[2:n:end];
ob3a = omg1[3:n:end];
ob4a = omg1[4:n:end];
ob5a = omg1[5:n:end];
ob6a = omg1[6:n:end];
ob7a = omg1[7:n:end];
ob8a = omg1[8:n:end];
ob9a = omg1[9:n:end];
ob10a = omg1[10:n:end];
ob11a = omg1[11:n:end];
ob12a = omg1[12:n:end];
ob13a = omg1[13:n:end];
ob14a = omg1[14:n:end];
ob15a = omg1[15:n:end];

gmxa = [maximum(gb1a[.!isnan.(gb1a)]),maximum(gb2a[.!isnan.(gb2a)]),maximum(gb3a[.!isnan.(gb3a)]),maximum(gb4a[.!isnan.(gb4a)]),maximum(gb5a[.!isnan.(gb5a)]),maximum(gb6a[.!isnan.(gb6a)]),maximum(gb7a[.!isnan.(gb7a)]),maximum(gb8a[.!isnan.(gb8a)]),maximum(gb9a[.!isnan.(gb9a)]),maximum(gb10a[.!isnan.(gb10a)]),maximum(gb11a[.!isnan.(gb11a)]),maximum(gb12a[.!isnan.(gb12a)]),maximum(gb13a[.!isnan.(gb13a)]),maximum(gb14a[.!isnan.(gb14a)]),maximum(gb15a[.!isnan.(gb15a)]),] 
# maxgam = convert(Array{Float64}, maxgam); 

function mxindxa(x)
 	for i in 1:28
 		if gmxa[x] == (gbxa[x])[i]
 		return i	
 		end
 	end
end

omxa = [ob1a[mxindxa(1)], ob2a[mxindxa(2)], ob3a[mxindxa(3)], ob4a[mxindxa(4)], ob5a[mxindxa(5)], ob6a[mxindxa(6)], ob7a[mxindxa(7)], ob8a[mxindxa(8)], ob9a[mxindxa(9)], ob10a[mxindxa(10)], ob11a[mxindxa(11)], ob12a[mxindxa(12)], ob13a[mxindxa(13)], ob14a[mxindxa(14)], ob15a[mxindxa(15)]]

gb1b = gam2[1:n:end];
gb2b = gam2[2:n:end];
gb3b = gam2[3:n:end];
gb4b = gam2[4:n:end];
gb5b = gam2[5:n:end];
gb6b = gam2[6:n:end];
gb7b = gam2[7:n:end];
gb8b = gam2[8:n:end];
gb9b = gam2[9:n:end];
gb10b = gam2[10:n:end];
gb11b = gam2[11:n:end];
gb12b = gam2[12:n:end];
gb13b = gam2[13:n:end];
gb14b = gam2[14:n:end];
gb15b = gam2[15:n:end];

gbxb = [gb1b, gb2b, gb3b, gb4b, gb5b, gb6b, gb7b, gb8b, gb9b, gb10b, gb11b, gb12b, gb13b, gb14b, gb15b]

ob1b = omg2[1:n:end];
ob2b = omg2[2:n:end];
ob3b = omg2[3:n:end];
ob4b = omg2[4:n:end];
ob5b = omg2[5:n:end];
ob6b = omg2[6:n:end];
ob7b = omg2[7:n:end];
ob8b = omg2[8:n:end];
ob9b = omg2[9:n:end];
ob10b = omg2[10:n:end];
ob11b = omg2[11:n:end];
ob12b = omg2[12:n:end];
ob13b = omg2[13:n:end];
ob14b = omg2[14:n:end];
ob15b = omg2[15:n:end];

gmxb = [maximum(gb1b[.!isnan.(gb1b)]),maximum(gb2b[.!isnan.(gb2b)]),maximum(gb3b[.!isnan.(gb3b)]),maximum(gb4b[.!isnan.(gb4b)]),maximum(gb5b[.!isnan.(gb5b)]),maximum(gb6b[.!isnan.(gb6b)]),maximum(gb7b[.!isnan.(gb7b)]),maximum(gb8b[.!isnan.(gb8b)]),maximum(gb9b[.!isnan.(gb9b)]),maximum(gb10b[.!isnan.(gb10b)]),maximum(gb11b[.!isnan.(gb11b)]),maximum(gb12b[.!isnan.(gb12b)]),maximum(gb13b[.!isnan.(gb13b)]),maximum(gb14b[.!isnan.(gb14b)]),maximum(gb15b[.!isnan.(gb15b)]),] 
# maxgam = convert(Array{Float64}, maxgam); 

function mxindxb(x)
 	for i in 1:28
 		if gmxb[x] == (gbxb[x])[i]
 		return i	
 		end
 	end
end

omxb = [ob1b[mxindxb(1)], ob2b[mxindxb(2)], ob3b[mxindxb(3)], ob4b[mxindxb(4)], ob5b[mxindxb(5)], ob6b[mxindxb(6)], ob7b[mxindxb(7)], ob8b[mxindxb(8)], ob9b[mxindxb(9)], ob10b[mxindxb(10)], ob11b[mxindxb(11)], ob12b[mxindxb(12)], ob13b[mxindxb(13)], ob14b[mxindxb(14)], ob15b[mxindxb(15)]]
# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"?? / %", ylabel = L"?? / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

gmx_itg = scatter!(beta7a, gmx[1:7], color = :green, markersize = markersize, marker = '???') 
gmx_kbm = scatter!(beta7b, gmx[8:15], color = :green, markersize = markersize, marker = '???') 
gmxl = lines!(ax1, beta7a, gmx[1:7], color = :green, linewidth = linewidth)
gmxl2 = lines!(ax1, beta7b, gmx[8:15], color = :green, linewidth = linewidth)

gmxa_itg = scatter!(beta8a, gmxa[1:8], color = :blue, markersize = markersize, marker = '???')  
gmxa_kbm = scatter!(beta8b, gmxa[9:15], color = :blue, markersize = markersize, marker = '???')
gmxal = lines!(ax1, beta8a, gmxa[1:8], color = :blue, linewidth = linewidth)
gmxal2 = lines!(ax1, beta8b, gmxa[9:15], color = :blue, linewidth = linewidth)

gmxb_itg = scatter!(beta6a, gmxb[1:6], color = :red, markersize = markersize, marker = '???')  
gmxb_kbm = scatter!(beta6b, gmxb[7:15], color = :red, markersize = markersize, marker = '???')
gmxbl = lines!(ax1, beta6a, gmxb[1:6], color = :red, linewidth = linewidth)
gmxbl2 = lines!(ax1, beta6b, gmxb[7:15], color = :red, linewidth = linewidth)

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth) 
hlines!(ax1, yax, color = :grey)
vlines!(ax1, xax, color = :grey)

ax2 = Axis(fig[1,2], xlabel = L"?? / %", ylabel = L"????? / (c_s/a)", xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

scatter!(beta7a, omx[1:7], color = :green, markersize = markersize, marker = '???')
scatter!(beta7b, omx[8:15], color = :green, markersize = markersize, marker = '???')
lines!(ax2, beta7a, omx[1:7], color = :green, linewidth = linewidth) 
lines!(ax2, beta7b, omx[8:15], color = :green, linewidth = linewidth)

scatter!(beta8a, omxa[1:8], color = :blue, markersize = markersize, marker = '???')
scatter!(beta8b, omxa[9:15], color = :blue, markersize = markersize, marker = '???')
lines!(ax2, beta8a, omxa[1:8], color = :blue, linewidth = linewidth) 
lines!(ax2, beta8b, omxa[9:15], color = :blue, linewidth = linewidth)

scatter!(beta6a, omxb[1:6], color = :red, markersize = markersize, marker = '???')
scatter!(beta6b, omxb[7:15], color = :red, markersize = markersize, marker = '???')
lines!(ax2, beta6a, omxb[1:6], color = :red, linewidth = linewidth) 
lines!(ax2, beta6b, omxb[7:15], color = :red, linewidth = linewidth)

hlines!(ax2, [0], color = :black, linewidth = linewidth) 
vlines!(ax2, [0], color = :black, linewidth = linewidth)
hlines!(ax2, yax2, color = :grey)
vlines!(ax2, xax, color = :grey)

xlims!(ax1, 0, 3.25) 
ylims!(ax1, 0, 0.45)
xlims!(ax2, 0, 3.25)
ylims!(ax2, 0, 0.65)


fig
