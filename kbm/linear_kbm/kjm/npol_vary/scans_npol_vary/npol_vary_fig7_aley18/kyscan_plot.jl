using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("scan.log_npol1");
sc1 = readdlm("scan.log_npol2");

n=8

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

gam1 = sc1[:,7];
gam1 = convert(Array{Float64}, gam1);

omg1 = sc1[:,8];
omg1 = convert(Array{Float64}, omg1);


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
gambet6 = gam[7:n:end];
gambet7 = gam[8:n:end];

omgbet0 = omg[1:n:end];
omgbet1 = omg[2:n:end];
omgbet2 = omg[3:n:end];
omgbet3 = omg[4:n:end];
omgbet4 = omg[5:n:end];
omgbet5 = omg[6:n:end];
omgbet6 = omg[6:n:end];
omgbet7 = omg[6:n:end];

gambet0a = gam1[1:n:end];
gambet1a = gam1[2:n:end];
gambet2a = gam1[3:n:end];
gambet3a = gam1[4:n:end];
gambet4a = gam1[5:n:end];
gambet5a = gam1[6:n:end];
gambet6a = gam1[7:n:end];
gambet7a = gam1[8:n:end];

omgbet0a = omg1[1:n:end];
omgbet1a = omg1[2:n:end];
omgbet2a = omg1[3:n:end];
omgbet3a = omg1[4:n:end];
omgbet4a = omg1[5:n:end];
omgbet5a = omg1[6:n:end];
omgbet6a = omg1[6:n:end];
omgbet7a = omg1[6:n:end];


# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])
ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

# hlines!(ax1, [gambet0[1], gambet1[1], gambet2[1], gambet3[1], gambet4[1], gambet5[1]], xmax = [1, 1, 1, 1, 1, 1], color = :black)
lines!(ax1, kyrho, gambet0, color = :red, label = "beta=0")
lines!(ax1, kyrho, gambet1, color = :orange, label = "beta=0.001")
lines!(ax1, kyrho, gambet2, color = :brown, label = "beta=0.002")
lines!(ax1, kyrho, gambet3, color = :green, label = "beta=0.003")
lines!(ax1, kyrho, gambet4, color = :cyan, label = "beta=0.004")
lines!(ax1, kyrho, gambet5, color = :blue, label = "beta=0.004")
lines!(ax1, kyrho, gambet6, color = :purple, label = "beta=0.004")
lines!(ax1, kyrho, gambet7, color = :black, label = "beta=0.004")

lines!(ax2, kyrho, omgbet0, color = :red, label = "beta=0")
lines!(ax2, kyrho, omgbet1, color = :orange, label = "beta=0.001")
lines!(ax2, kyrho, omgbet2, color = :brown, label = "beta=0.002")
lines!(ax2, kyrho, omgbet3, color = :green, label = "beta=0.003")
lines!(ax2, kyrho, omgbet4, color = :cyan, label = "beta=0.004")
lines!(ax2, kyrho, omgbet5, color = :blue, label = "beta=0.005")
lines!(ax2, kyrho, omgbet6, color = :purple, label = "beta=0.005")
lines!(ax2, kyrho, omgbet7, color = :black, label = "beta=0.005")

lines!(ax3, kyrho, gambet0a, color = :red, label = "beta=0")
lines!(ax3, kyrho, gambet1a, color = :orange, label = "beta=0.001")
lines!(ax3, kyrho, gambet2a, color = :brown, label = "beta=0.002")
lines!(ax3, kyrho, gambet3a, color = :green, label = "beta=0.003")
lines!(ax3, kyrho, gambet4a, color = :cyan, label = "beta=0.004")
lines!(ax3, kyrho, gambet5a, color = :blue, label = "beta=0.004")
lines!(ax3, kyrho, gambet6a, color = :purple, label = "beta=0.004")
lines!(ax3, kyrho, gambet7a, color = :black, label = "beta=0.004")

lines!(ax4, kyrho, omgbet0a, color = :red, label = "beta=0")
lines!(ax4, kyrho, omgbet1a, color = :orange, label = "beta=0.001")
lines!(ax4, kyrho, omgbet2a, color = :brown, label = "beta=0.002")
lines!(ax4, kyrho, omgbet3a, color = :green, label = "beta=0.003")
lines!(ax4, kyrho, omgbet4a, color = :cyan, label = "beta=0.004")
lines!(ax4, kyrho, omgbet5a, color = :blue, label = "beta=0.005")
lines!(ax4, kyrho, omgbet6a, color = :purple, label = "beta=0.005")
lines!(ax4, kyrho, omgbet7a, color = :black, label = "beta=0.005")

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

ax3.title = "γ"
ax4.title = "ωᵣ"

ax3.xlabel = "kyρ"
ax4.xlabel = "kyρ"


#xlims!(ax1, 0, 0.90)
#xlims!(ax2, 0, 0.90)
#ylims!(ax1, 0, 0.70)
#ylims!(ax2, 0, 0.80)

save("evalue_kyscan_npol1.png",fig)

fig
