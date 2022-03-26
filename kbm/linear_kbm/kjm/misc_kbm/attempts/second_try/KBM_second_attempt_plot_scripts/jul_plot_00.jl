using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("00scan");

beta = sc[1:6,3];
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

gambet0 = gam[1:6:end];
gambet1 = gam[2:6:end];
gambet2 = gam[3:6:end];
gambet3 = gam[4:6:end];
gambet4 = gam[5:6:end];
gambet5 = gam[6:6:end];

omgbet0 = omg[1:6:end];
omgbet1 = omg[2:6:end];
omgbet2 = omg[3:6:end];
omgbet3 = omg[4:6:end];
omgbet4 = omg[5:6:end];
omgbet5 = omg[6:6:end];
# plot(kyrho, gambet0)
fig = Figure()

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

lines!(ax1, kyrho, gambet0, color = :red, label = "beta=0")
lines!(ax1, kyrho, gambet1, color = :orange, label = "beta=0.001")
lines!(ax1, kyrho, gambet2, color = :yellow, label = "beta=0.002")
lines!(ax1, kyrho, gambet3, color = :green, label = "beta=0.003")
lines!(ax1, kyrho, gambet4, color = :blue, label = "beta=0.004")
lines!(ax1, kyrho, gambet5, color = :purple, label = "beta=0.005")

lines!(ax2, kyrho, omgbet0, color = :red, label = "beta=0")
lines!(ax2, kyrho, omgbet1, color = :orange, label = "beta=0.001")
lines!(ax2, kyrho, omgbet2, color = :yellow, label = "beta=0.002")
lines!(ax2, kyrho, omgbet3, color = :green, label = "beta=0.003")
lines!(ax2, kyrho, omgbet4, color = :blue, label = "beta=0.004")
lines!(ax2, kyrho, omgbet5, color = :purple, label = "beta=0.005")

ax1.title = "Gamma"
ax2.title = "Omega_r"

ax1.xlabel = "kyrho"
ax2.xlabel = "kyrho"

save("evalue_scan_00.png",fig)

fig
