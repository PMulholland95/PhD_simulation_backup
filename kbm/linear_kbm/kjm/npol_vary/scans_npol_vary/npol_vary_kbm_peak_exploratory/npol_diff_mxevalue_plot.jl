using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("edit_npol1")
sc1 = readdlm("edit_npol2")
n=11

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

gam1 = sc1[:,7];
gam1 = convert(Array{Float64}, gam1);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

omg1 = sc1[:,8];
omg1 = convert(Array{Float64}, omg1);

kyrho= [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09];

# We now want to gather together gamma values for range of kyrho, and fixed beta
# In the end, only want to plot: gamma|max vs. beta 

#Makes array (gambet0) of every 6th value in original array (gam) starting from first value (beta=0)
# This gathers all gammas with beta=0 for all kyrho values in scan

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

gbx = [gb1, gb2, gb3, gb4, gb5, gb6, gb7, gb8, gb9, gb10, gb11]

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

gbxa = [gb1a, gb2a, gb3a, gb4a, gb5a, gb6a, gb7a, gb8a, gb9a, gb10a, gb11a]

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

gmx = [maximum(gb1[.!isnan.(gb1)]),maximum(gb2[.!isnan.(gb2)]),maximum(gb3[.!isnan.(gb3)]),maximum(gb4[.!isnan.(gb4)]),maximum(gb5[.!isnan.(gb5)]),maximum(gb6[.!isnan.(gb6)]),maximum(gb7[.!isnan.(gb7)]),maximum(gb8[.!isnan.(gb8)]),maximum(gb9[.!isnan.(gb9)]),maximum(gb10[.!isnan.(gb10)]),maximum(gb11[.!isnan.(gb11)]),] 
# maxgam = convert(Array{Float64}, maxgam); 

gmxa = [maximum(gb1a[.!isnan.(gb1a)]),maximum(gb2a[.!isnan.(gb2a)]),maximum(gb3a[.!isnan.(gb3a)]),maximum(gb4a[.!isnan.(gb4a)]),maximum(gb5a[.!isnan.(gb5a)]),maximum(gb6a[.!isnan.(gb6a)]),maximum(gb7a[.!isnan.(gb7a)]),maximum(gb8a[.!isnan.(gb8a)]),maximum(gb9a[.!isnan.(gb9a)]),maximum(gb10a[.!isnan.(gb10a)]),maximum(gb11a[.!isnan.(gb11a)]),] 

gdif0 = gmx - gmxa
gdif = map(x->abs((100*gdif0[x])/(gmx[x])), 1:length(gdif0))

function mxindx(x)
 	for i in 1:10
 		if gmx[x] == (gbx[x])[i]
 		return i	
 		end
 	end
end

omx = [ob1[mxindx(1)], ob2[mxindx(2)], ob3[mxindx(3)], ob4[mxindx(4)], ob5[mxindx(5)], ob6[mxindx(6)], ob7[mxindx(7)], ob8[mxindx(8)], ob9[mxindx(9)], ob10[mxindx(10)], ob11[mxindx(11)]]

function mxindxa(x)
 	for i in 1:10
 		if gmxa[x] == (gbxa[x])[i]
 		return i	
 		end
 	end
end

omxa = [ob1a[mxindxa(1)], ob2a[mxindxa(2)], ob3a[mxindxa(3)], ob4a[mxindxa(4)], ob5a[mxindxa(5)], ob6a[mxindxa(6)], ob7a[mxindxa(7)], ob8a[mxindxa(8)], ob9a[mxindxa(9)], ob10a[mxindxa(10)], ob11a[mxindxa(11)]]

odif0 = omx - omxa
odif = map(x->abs((100*odif0[x])/(omx[x])), 1:length(odif0))

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

gmxl = lines!(ax1, beta, gmx, color = :red, label = "beta=0")
gmxal = lines!(ax1, beta, gmxa, color = :green, label = "beta=0")

omxl = lines!(ax2, beta, omx, color = :blue, label = "beta=0")
omxal = lines!(ax2, beta, omxa, color = :orange, label = "beta=0")

ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

lines!(ax3, beta, gmxa, color = :green, label = "beta=0")

lines!(ax4, beta, omxa, color = :orange, label = "beta=0")

ax5 = Axis(fig[3,1])
ax6 = Axis(fig[3,2])

lines!(ax5, beta, gdif, color = :black, label = "beta=0")

lines!(ax6, beta, odif, color = :black, label = "beta=0")

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "β"
ax2.xlabel = "β"

ax3.title = "γ"
ax4.title = "ωᵣ"

ax3.xlabel = "β"
ax4.xlabel = "β"

ax5.title = "γ (% difference b/t npol=1,2)"
ax6.title = "ωᵣ (% difference b/t npol=1,2)"

ax5.xlabel = "β"
ax6.xlabel = "β"

Legend(fig[1,3], [gmxl, omxl, gmxal, omxal], ["npol = 1", "npol = 1", "npol = 2", "npol = 2"])
#xlims!(ax1, 0.0, 0.037)
#xlims!(ax2, 0.0, 0.037)
#ylims!(ax1, 0, 0.5)
#ylims!(ax2, 0, 0.5)

save("npol_diff_maxevalue_betascan.png",fig)

fig
