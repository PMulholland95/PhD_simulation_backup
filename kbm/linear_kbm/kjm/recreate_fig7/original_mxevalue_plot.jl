using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("scan.log");

n=8

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

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

gmx = [maximum(gb0[.!isnan.(gb0)]),maximum(gb1[.!isnan.(gb1)]),maximum(gb2[.!isnan.(gb2)]),maximum(gb3[.!isnan.(gb3)]),maximum(gb4[.!isnan.(gb4)]),maximum(gb5[.!isnan.(gb5)]),maximum(gb6[.!isnan.(gb6)]),maximum(gb7[.!isnan.(gb7)])] 
# maxgam = convert(Array{Float64}, maxgam); 

omx = [ob0[16], ob1[16], ob2[16], ob3[16], ob4[16], ob5[1], ob6[1], ob7[1]]

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

lines!(ax1, beta, gmx, color = :red, label = "beta=0")

lines!(ax2, beta, omx, color = :blue, label = "beta=0")

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "β"
ax2.xlabel = "β"

xlims!(ax1, 0.0, 0.037)
xlims!(ax2, 0.0, 0.037)
ylims!(ax1, 0, 0.5)
ylims!(ax2, 0, 0.5)

save("maxevalue_betascan_45.png",fig)

fig
