using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

#sc = readdlm("omtn_0_coarse_scan");
sc0 = readdlm("omtn_0_coarse_scan")
sc1 = readdlm("omt_1_omn_0_coarse_scan")
sc2 = readdlm("omt_0_omn_1_coarse_scan")
sc3 = readdlm("omt_1_omn_1_coarse_scan")

n=8

beta = sc0[1:n,3];
beta = convert(Array{Float64}, beta);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);

gam1 = sc1[:,7];
gam1 = convert(Array{Float64}, gam1);
omg1 = sc1[:,8];
omg1 = convert(Array{Float64}, omg1);

gam2 = sc2[:,7];
gam2 = convert(Array{Float64}, gam2);
omg2 = sc2[:,8];
omg2 = convert(Array{Float64}, omg2);

gam3 = sc3[:,7];
gam3 = convert(Array{Float64}, gam3);
omg3 = sc3[:,8];
omg3 = convert(Array{Float64}, omg3);

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

gb0b = gam1[1:n:end];
gb1b = gam1[2:n:end];
gb2b = gam1[3:n:end];
gb3b = gam1[4:n:end];
gb4b = gam1[5:n:end];
gb5b = gam1[6:n:end];
gb6b = gam1[7:n:end];
gb7b = gam1[8:n:end];

ob0b = omg1[1:n:end];
ob1b = omg1[2:n:end];
ob2b = omg1[3:n:end];
ob3b = omg1[4:n:end];
ob4b = omg1[5:n:end];
ob5b = omg1[6:n:end];
ob6b = omg1[7:n:end];
ob7b = omg1[8:n:end];

gb0c = gam2[1:n:end];
gb1c = gam2[2:n:end];
gb2c = gam2[3:n:end];
gb3c = gam2[4:n:end];
gb4c = gam2[5:n:end];
gb5c = gam2[6:n:end];
gb6c = gam2[7:n:end];
gb7c = gam2[8:n:end];

ob0c = omg2[1:n:end];
ob1c = omg2[2:n:end];
ob2c = omg2[3:n:end];
ob3c = omg2[4:n:end];
ob4c = omg2[5:n:end];
ob5c = omg2[6:n:end];
ob6c = omg2[7:n:end];
ob7c = omg2[8:n:end];

gb0d = gam3[1:n:end];
gb1d = gam3[2:n:end];
gb2d = gam3[3:n:end];
gb3d = gam3[4:n:end];
gb4d = gam3[5:n:end];
gb5d = gam3[6:n:end];
gb6d = gam3[7:n:end];
gb7d = gam3[8:n:end];

ob0d = omg3[1:n:end];
ob1d = omg3[2:n:end];
ob2d = omg3[3:n:end];
ob3d = omg3[4:n:end];
ob4d = omg3[5:n:end];
ob5d = omg3[6:n:end];
ob6d = omg3[7:n:end];
ob7d = omg3[8:n:end];

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

ax3 = Axis(fig[2,1])
ax4 = Axis(fig[2,2])

ax5 = Axis(fig[3,1])
ax6 = Axis(fig[3,2])

ax7 = Axis(fig[4,1])
ax8 = Axis(fig[4,2])

#hlines!(ax1, [gambet0[1], gambet1[1], gambet2[1], gambet3[1], gambet4[1], gambet5[1]], xmax = [1, 1, 1, 1, 1, 1], color = :black)
gblines0 = lines!(ax1, kyrho, gb0a, color = :red)
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

lines!(ax3, kyrho, gb0b, color = :red)
lines!(ax3, kyrho, gb1b, color = :orange)
lines!(ax3, kyrho, gb2b, color = :brown)
lines!(ax3, kyrho, gb3b, color = :green)
lines!(ax3, kyrho, gb4b, color = :cyan)
lines!(ax3, kyrho, gb5b, color = :blue)
lines!(ax3, kyrho, gb6b, color = :purple)
lines!(ax3, kyrho, gb7b, color = :black) 

lines!(ax4, kyrho, ob0b, color = :red)
lines!(ax4, kyrho, ob1b, color = :orange)
lines!(ax4, kyrho, ob2b, color = :brown)
lines!(ax4, kyrho, ob3b, color = :green)
lines!(ax4, kyrho, ob4b, color = :cyan)
lines!(ax4, kyrho, ob5b, color = :blue)
lines!(ax4, kyrho, ob6b, color = :purple)
lines!(ax4, kyrho, ob7b, color = :black) 

lines!(ax5, kyrho, gb0c, color = :red)
lines!(ax5, kyrho, gb1c, color = :orange)
lines!(ax5, kyrho, gb2c, color = :brown)
lines!(ax5, kyrho, gb3c, color = :green)
lines!(ax5, kyrho, gb4c, color = :cyan)
lines!(ax5, kyrho, gb5c, color = :blue)
lines!(ax5, kyrho, gb6c, color = :purple)
lines!(ax5, kyrho, gb7c, color = :black) 

lines!(ax6, kyrho, ob0c, color = :red)
lines!(ax6, kyrho, ob1c, color = :orange)
lines!(ax6, kyrho, ob2c, color = :brown)
lines!(ax6, kyrho, ob3c, color = :green)
lines!(ax6, kyrho, ob4c, color = :cyan)
lines!(ax6, kyrho, ob5c, color = :blue)
lines!(ax6, kyrho, ob6c, color = :purple)

lines!(ax7, kyrho, gb0d, color = :red)
lines!(ax7, kyrho, gb1d, color = :orange)
lines!(ax7, kyrho, gb2d, color = :brown)
lines!(ax7, kyrho, gb3d, color = :green)
lines!(ax7, kyrho, gb4d, color = :cyan)
lines!(ax7, kyrho, gb5d, color = :blue)
lines!(ax7, kyrho, gb6d, color = :purple)
lines!(ax7, kyrho, gb7d, color = :black) 

lines!(ax8, kyrho, ob0d, color = :red)
lines!(ax8, kyrho, ob1d, color = :orange)
lines!(ax8, kyrho, ob2d, color = :brown)
lines!(ax8, kyrho, ob3d, color = :green)
lines!(ax8, kyrho, ob4d, color = :cyan)
lines!(ax8, kyrho, ob5d, color = :blue)
lines!(ax8, kyrho, ob6d, color = :purple)

ax1.title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=0 ; omnₑ=0)"
ax2.title = "ωᵣ"

ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

ax3.title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=1 ; omnₑ=0)"
ax4.title = "ωᵣ"

ax3.xlabel = "kyρ"
ax4.xlabel = "kyρ"

ax5.title = "γ (omtᵢ=3 ; omnᵢ=1 ; omtₑ=0 ; omnₑ=1)"
ax6.title = "ωᵣ"

ax5.xlabel = "kyρ"
ax6.xlabel = "kyρ"

ax7.title = "γ (omtᵢ=3 ; omnᵢ=1 ; omtₑ=1 ; omnₑ=1)"
ax8.title = "ωᵣ"

ax7.xlabel = "kyρ"
ax8.xlabel = "kyρ"

Legend(fig[1, 3],
	[gblines0, gblines1, gblines2, gblines3, gblines4, gblines5, gblines6, gblines7],
	["β = 0.0%","β = 0.5%","β = 1.0%","β = 1.5%","β = 2.0%","β = 2.5%","β = 3.0%","β = 3.5%"])


#xlims!(ax1, 0, 0.16)
#xlims!(ax2, 0, 0.16)
#ylims!(ax1, 0, 0.45)
#ylims!(ax2, 0, 0.32)

save("evalue_coarse_kyscan.png",fig)

fig
