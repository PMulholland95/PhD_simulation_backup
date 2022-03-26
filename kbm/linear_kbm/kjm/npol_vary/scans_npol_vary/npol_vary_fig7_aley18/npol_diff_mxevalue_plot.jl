using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("scan.log_npol1")
sc1 = readdlm("scan.log_npol2")
n=8

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

kyrho=range(0.05,0.8, length=16);

# We now want to gather together gamma values for range of kyrho, and fixed beta
# In the end, only want to plot: gamma|max vs. beta 

#Makes array (gambet0) of every 6th value in original array (gam) starting from first value (beta=0)
# This gathers all gammas with beta=0 for all kyrho values in scan

gb0 = gam[1:n:end];
gb1 = gam[2:n:end];
gb2 = gam[3:n:end];
gb3 = gam[4:n:end];
gb4 = gam[5:n:end];
gb5 = gam[6:n:end];
gb6 = gam[7:n:end];
gb7 = gam[8:n:end];

ob0 = omg[1:n:end];
ob1 = omg[2:n:end];
ob2 = omg[3:n:end];
ob3 = omg[4:n:end];
ob4 = omg[5:n:end];
ob5 = omg[6:n:end];
ob6 = omg[6:n:end];
ob7 = omg[6:n:end];

gb0a = gam1[1:n:end];
gb1a = gam1[2:n:end];
gb2a = gam1[3:n:end];
gb3a = gam1[4:n:end];
gb4a = gam1[5:n:end];
gb5a = gam1[6:n:end];
gb6a = gam1[7:n:end];
gb7a = gam1[8:n:end];

ob0a = omg1[1:n:end];
ob1a = omg1[2:n:end];
ob2a = omg1[3:n:end];
ob3a = omg1[4:n:end];
ob4a = omg1[5:n:end];
ob5a = omg1[6:n:end];
ob6a = omg1[6:n:end];
ob7a = omg1[6:n:end];

gmx = [maximum(gb0[.!isnan.(gb0)]),maximum(gb1[.!isnan.(gb1)]),maximum(gb2[.!isnan.(gb2)]),maximum(gb3[.!isnan.(gb3)]),maximum(gb4[.!isnan.(gb4)]),maximum(gb5[.!isnan.(gb5)]),maximum(gb6[.!isnan.(gb6)]),maximum(gb7[.!isnan.(gb7)])] 
# maxgam = convert(Array{Float64}, maxgam); 

gmxa = [maximum(gb0a[.!isnan.(gb0a)]),maximum(gb1a[.!isnan.(gb1a)]),maximum(gb2a[.!isnan.(gb2a)]),maximum(gb3a[.!isnan.(gb3a)]),maximum(gb4a[.!isnan.(gb4a)]),maximum(gb5a[.!isnan.(gb5a)]),maximum(gb6a[.!isnan.(gb6a)]),maximum(gb7a[.!isnan.(gb7a)])] 

gdif0 = gmx - gmxa
gdif = map(x->abs((100*gdif0[x])/(gmx[x])), 1:length(gdif0))

omx = [ob0[16], ob1[16], ob2[1], ob3[1], ob4[1], ob5[1], ob6[1], ob7[1]]

omxa = [ob0a[16], ob1a[16], ob2a[1], ob3a[1], ob4a[1], ob5a[1], ob6a[1], ob7a[1]]

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

save("maxevalue_betascan.png",fig)

fig
