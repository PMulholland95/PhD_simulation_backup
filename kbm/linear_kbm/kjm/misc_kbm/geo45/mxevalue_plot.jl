using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("45scan");

n=6

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

gambet0 = gam[1:n:end];
gambet1 = gam[2:n:end];
gambet2 = gam[3:n:end];
gambet3 = gam[4:n:end];
gambet4 = gam[5:n:end];
gambet5 = gam[6:n:end];

omgbet0 = omg[1:n:end];
omgbet1 = omg[2:n:end];
omgbet2 = omg[3:n:end];
omgbet3 = omg[4:n:end];
omgbet4 = omg[5:n:end];
omgbet5 = omg[6:n:end];

gmx = [maximum(gambet0),maximum(gambet1),maximum(gambet2),maximum(gambet3),maximum(gambet4),maximum(gambet5)] 
# maxgam = convert(Array{Float64}, maxgam); 

omx = [maximum(omgbet0),maximum(omgbet1),maximum(omgbet2),maximum(omgbet3),maximum(omgbet4),maximum(omgbet5)] 
# maxomg = convert(Array{Float64}, maxomg); 

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

xlims!(ax1, 0.024, 0.029)
xlims!(ax2, 0.024, 0.029)
ylims!(ax1, 0, 0.25)
ylims!(ax2, 0, 1.00)

save("maxevalue_betascan_45.png",fig)

fig
