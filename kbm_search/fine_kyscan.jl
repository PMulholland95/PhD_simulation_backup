using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

#sc = readdlm("omtn_0_fine_scan");
sc = readdlm("omt_1_omn_0_fine_scan");

n=11

beta = sc[1:n,3];
beta = convert(Array{Float64}, beta);

gam = sc[:,7];
gam = convert(Array{Float64}, gam);

omg = sc[:,8];
omg = convert(Array{Float64}, omg);

# kyrho=range(1,5, length=5);
kyrho= [0.01, 0.02, 0.03, 0.04, 0.05]

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
gambet8 = gam[9:n:end];
gambet9 = gam[10:n:end];
gambet10 = gam[11:n:end];

omgbet0 = omg[1:n:end];
omgbet1 = omg[2:n:end];
omgbet2 = omg[3:n:end];
omgbet3 = omg[4:n:end];
omgbet4 = omg[5:n:end];
omgbet5 = omg[6:n:end];
omgbet6 = omg[7:n:end];
omgbet7 = omg[8:n:end];
omgbet8 = omg[9:n:end];
omgbet9 = omg[10:n:end];
omgbet10 = omg[11:n:end];

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 30) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

#hlines!(ax1, [gambet0[1], gambet1[1], gambet2[1], gambet3[1], gambet4[1], gambet5[1]], xmax = [1, 1, 1, 1, 1, 1], color = :black)
gblines0 = lines!(ax1, kyrho, gambet0, color = :red)
gblines1 = lines!(ax1, kyrho, gambet1, color = :orange)
gblines2 = lines!(ax1, kyrho, gambet2, color = :brown)
gblines3 = lines!(ax1, kyrho, gambet3, color = :green)
gblines4 = lines!(ax1, kyrho, gambet4, color = :cyan)
gblines5 = lines!(ax1, kyrho, gambet5, color = :blue)
gblines6 = lines!(ax1, kyrho, gambet6, color = :purple)
gblines7 = lines!(ax1, kyrho, gambet7, color = :black)
gblines8 = lines!(ax1, kyrho, gambet8, color = :red)
gblines9 = lines!(ax1, kyrho, gambet9, color = :orange)
gblines10 = lines!(ax1, kyrho, gambet10, color = :brown)

lines!(ax2, kyrho, omgbet0, color = :red)
lines!(ax2, kyrho, omgbet1, color = :orange)
lines!(ax2, kyrho, omgbet2, color = :brown)
lines!(ax2, kyrho, omgbet3, color = :green)
lines!(ax2, kyrho, omgbet4, color = :cyan)
lines!(ax2, kyrho, omgbet5, color = :blue)
lines!(ax2, kyrho, omgbet6, color = :purple)
lines!(ax2, kyrho, omgbet7, color = :black)
lines!(ax2, kyrho, omgbet8, color = :red)
lines!(ax2, kyrho, omgbet9, color = :orange)
lines!(ax2, kyrho, omgbet10, color = :brown)

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

Legend(fig[1, 3],
       [gblines0, gblines1, gblines2, gblines3, gblines4, gblines5, gblines6, gblines7, gblines8, gblines9, gblines10],
       ["β = 0.0%","β = 0.2%","β = 0.4%","β = 0.6%","β = 0.8%","β = 1.0%","β = 1.2%","β = 1.4%","β = 1.6%","β = 1.8%","β = 2.0%",])

xlims!(ax1, 0, 0.06)
xlims!(ax2, 0, 0.06)
ylims!(ax1, 0, 0.30)
ylims!(ax2, 0, 0.10)

save("evalue_fine_kyscan.png",fig)

fig
