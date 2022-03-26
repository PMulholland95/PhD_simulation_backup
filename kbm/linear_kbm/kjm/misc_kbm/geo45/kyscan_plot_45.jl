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

maxgam = (maximum(gambet0),maximum(gambet1),maximum(gambet2),maximum(gambet3),maximum(gambet4),maximum(gambet5)) 

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

hlines!(ax1, [gambet0[1], gambet1[1], gambet2[1], gambet3[1], gambet4[1], gambet5[1]], xmax = [1, 1, 1, 1, 1, 1], color = :black)
lines!(ax1, kyrho, gambet0, color = :red, label = "beta=0")
lines!(ax1, kyrho, gambet1, color = :orange, label = "beta=0.001")
lines!(ax1, kyrho, gambet2, color = :yellow, label = "beta=0.002")
lines!(ax1, kyrho, gambet3, color = :green, label = "beta=0.003")
lines!(ax1, kyrho, gambet4, color = :cyan, label = "beta=0.004")
lines!(ax1, kyrho, gambet5, color = :blue, label = "beta=0.004")

lines!(ax2, kyrho, omgbet0, color = :red, label = "beta=0")
lines!(ax2, kyrho, omgbet1, color = :orange, label = "beta=0.001")
lines!(ax2, kyrho, omgbet2, color = :yellow, label = "beta=0.002")
lines!(ax2, kyrho, omgbet3, color = :green, label = "beta=0.003")
lines!(ax2, kyrho, omgbet4, color = :cyan, label = "beta=0.004")
lines!(ax2, kyrho, omgbet5, color = :blue, label = "beta=0.005")

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

xlims!(ax1, 0, 0.90)
xlims!(ax2, 0, 0.90)
ylims!(ax1, 0, 0.25)
ylims!(ax2, 0, 1.00)

save("evalue_kyscan_45.png",fig)

fig
