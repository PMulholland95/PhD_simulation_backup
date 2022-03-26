using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

#sc = readdlm("omtn_0_coarse_scan");
sc0 = readdlm("coarse_gista_scan")

n=8

beta = sc0[1:n,3];
beta = convert(Array{Float64}, beta);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);


# kylist = 0.010, 0.025, 0.050, 0.100, 0.150
# kyrho=range(1,5, length=5);
kyrho= [0.010, 0.025, 0.050, 0.100, 0.150]

# We now want to gather together gamma values for range of kyrho, and fixed beta
# In the end, only want to plot: gamma|max vs. beta 

#Makes array (gambet0) of every 6th value in original array (gam) starting from first value (beta=0)
# This gathers all gammas with beta=0 for all kyrho values in scan

gb0a = gam0[1:n:end];
gb1a = gam0[2:n:end];
gb2a = gam0[3:n:end];
gb3a = gam0[4:n:end];
gb4a = gam0[5:n:end];
gb5a = gam0[6:n:end];
gb6a = gam0[7:n:end];
gb7a = gam0[8:n:end];

ob0a = omg0[1:n:end];
ob1a = omg0[2:n:end];
ob2a = omg0[3:n:end];
ob3a = omg0[4:n:end];
ob4a = omg0[5:n:end];
ob5a = omg0[6:n:end];
ob6a = omg0[7:n:end];
ob7a = omg0[8:n:end];

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

gblines1 = lines!(ax1, kyrho, gb1a, color = :orange)
gblines2 = lines!(ax1, kyrho, gb2a, color = :brown)
gblines3 = lines!(ax1, kyrho, gb3a, color = :green)
gblines4 = lines!(ax1, kyrho, gb4a, color = :cyan)
gblines5 = lines!(ax1, kyrho, gb5a, color = :blue)
gblines6 = lines!(ax1, kyrho, gb6a, color = :purple)
gblines7 = lines!(ax1, kyrho, gb7a, color = :black)

lines!(ax2, kyrho, ob0a, color = :red)
lines!(ax2, kyrho, ob1a, color = :orange)
lines!(ax2, kyrho, ob2a, color = :brown)
lines!(ax2, kyrho, ob3a, color = :green)
lines!(ax2, kyrho, ob4a, color = :cyan)
lines!(ax2, kyrho, ob5a, color = :blue)
lines!(ax2, kyrho, ob6a, color = :purple)
lines!(ax2, kyrho, ob7a, color = :black) 


ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

ax1.title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=0 ; omnₑ=0)"
ax2.title = "ωᵣ"

Legend(fig[1, 3],
	[gblines0, gblines1, gblines2, gblines3, gblines4, gblines5, gblines6, gblines7],
	["β = 0.0%","β = 0.5%","β = 1.0%","β = 1.5%","β = 2.0%","β = 2.5%","β = 3.0%","β = 3.5%"])


#xlims!(ax1, 0, 0.16)
#xlims!(ax2, 0, 0.16)
#ylims!(ax1, 0, 0.45)
#ylims!(ax2, 0, 0.32)

save("evalue_coarse_kyscan.png",fig)

fig
