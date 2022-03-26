using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

#sc = readdlm("omtn_0_fine_scan");
sc1 = readdlm("npol_1_omt_i_1.6_p16sc1");
sc2 = readdlm("npol_2_omt_i_1.6_p16sc2");
sc3 = readdlm("npol_3_nx0_60_omt_i_1.6_p16sc4");

n=11

beta = sc1[1:n,3];
beta = convert(Array{Float64}, beta);

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

# kyrho=range(1,5, length=5);
kyrho= [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.5]

gambet1 = gam1[1:n:end];
gambet2 = gam1[2:n:end];
gambet3 = gam1[3:n:end];
gambet4 = gam1[4:n:end];
gambet5 = gam1[5:n:end];
gambet6 = gam1[6:n:end];
gambet7 = gam1[7:n:end];
gambet8 = gam1[8:n:end];
gambet9 = gam1[9:n:end];
gambet10 = gam1[10:n:end];
gambet11 = gam1[11:n:end];

omgbet1 = omg1[1:n:end];
omgbet2 = omg1[2:n:end];
omgbet3 = omg1[3:n:end];
omgbet4 = omg1[4:n:end];
omgbet5 = omg1[5:n:end];
omgbet6 = omg1[6:n:end];
omgbet7 = omg1[7:n:end];
omgbet8 = omg1[8:n:end];
omgbet9 = omg1[9:n:end];
omgbet10 = omg1[10:n:end];
omgbet11 = omg1[11:n:end];

gambetA1 = gam2[1:n:end];
gambetA2 = gam2[2:n:end];
gambetA3 = gam2[3:n:end];
gambetA4 = gam2[4:n:end];
gambetA5 = gam2[5:n:end];
gambetA6 = gam2[6:n:end];
gambetA7 = gam2[7:n:end];
gambetA8 = gam2[8:n:end];
gambetA9 = gam2[9:n:end];
gambetA10 = gam2[10:n:end];
gambetA11 = gam2[11:n:end];

omgbetA1 = omg2[1:n:end];
omgbetA2 = omg2[2:n:end];
omgbetA3 = omg2[3:n:end];
omgbetA4 = omg2[4:n:end];
omgbetA5 = omg2[5:n:end];
omgbetA6 = omg2[6:n:end];
omgbetA7 = omg2[7:n:end];
omgbetA8 = omg2[8:n:end];
omgbetA9 = omg2[9:n:end];
omgbetA10 = omg2[10:n:end];
omgbetA11 = omg2[11:n:end];

gambetB1 = gam3[1:n:end];
gambetB2 = gam3[2:n:end];
gambetB3 = gam3[3:n:end];
gambetB4 = gam3[4:n:end];
gambetB5 = gam3[5:n:end];
gambetB6 = gam3[6:n:end];
gambetB7 = gam3[7:n:end];
gambetB8 = gam3[8:n:end];
gambetB9 = gam3[9:n:end];
gambetB10 = gam3[10:n:end];
gambetB11 = gam3[11:n:end];

omgbetB1 = omg3[1:n:end];
omgbetB2 = omg3[2:n:end];
omgbetB3 = omg3[3:n:end];
omgbetB4 = omg3[4:n:end];
omgbetB5 = omg3[5:n:end];
omgbetB6 = omg3[6:n:end];
omgbetB7 = omg3[7:n:end];
omgbetB8 = omg3[8:n:end];
omgbetB9 = omg3[9:n:end];
omgbetB10 = omg3[10:n:end];
omgbetB11 = omg3[11:n:end];

linewidth = 4 

# plot(kyrho, gambet0)
fig = Figure()

fontsize_theme = Theme(fontsize = 60) 
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

#hlines!(ax1, [gambet0[1], gambet1[1], gambet2[1], gambet3[1], gambet4[1], gambet5[1]], xmax = [1, 1, 1, 1, 1, 1], color = :black)
gblines4 = lines!(ax1, kyrho, gambet4,linewidth = linewidth, color = :cyan)
gblines5 = lines!(ax1, kyrho, gambet5,linewidth = linewidth, color = :blue)
gblines6 = lines!(ax1, kyrho, gambet6,linewidth = linewidth, color = :purple)
gblines7 = lines!(ax1, kyrho, gambet7,linewidth = linewidth, color = :black)

lines!(ax2, kyrho, omgbet4,linewidth = linewidth, color = :cyan)
lines!(ax2, kyrho, omgbet5,linewidth = linewidth, color = :blue)
lines!(ax2, kyrho, omgbet6,linewidth = linewidth, color = :purple)
lines!(ax2, kyrho, omgbet7,linewidth = linewidth, color = :black)

lines!(ax1, kyrho, gambetA4,linewidth = linewidth, color = :cyan, linestyle = :dash)
lines!(ax1, kyrho, gambetA5,linewidth = linewidth, color = :blue, linestyle = :dash)
lines!(ax1, kyrho, gambetA6,linewidth = linewidth, color = :purple, linestyle = :dash)
lines!(ax1, kyrho, gambetA7,linewidth = linewidth, color = :black, linestyle = :dash)

lines!(ax2, kyrho, omgbetA4,linewidth = linewidth, color = :cyan, linestyle = :dash)
lines!(ax2, kyrho, omgbetA5,linewidth = linewidth, color = :blue, linestyle = :dash)
lines!(ax2, kyrho, omgbetA6,linewidth = linewidth, color = :purple, linestyle = :dash)
lines!(ax2, kyrho, omgbetA7,linewidth = linewidth, color = :black, linestyle = :dash)

lines!(ax1, kyrho, gambetB4,linewidth = linewidth, color = :cyan, linestyle = :dot)
lines!(ax1, kyrho, gambetB5,linewidth = linewidth, color = :blue, linestyle = :dot)
lines!(ax1, kyrho, gambetB6,linewidth = linewidth, color = :purple, linestyle = :dot)
lines!(ax1, kyrho, gambetB7,linewidth = linewidth, color = :black, linestyle = :dot)

lines!(ax2, kyrho, omgbetB4,linewidth = linewidth, color = :cyan, linestyle = :dot)
lines!(ax2, kyrho, omgbetB5,linewidth = linewidth, color = :blue, linestyle = :dot)
lines!(ax2, kyrho, omgbetB6,linewidth = linewidth, color = :purple, linestyle = :dot)
lines!(ax2, kyrho, omgbetB7,linewidth = linewidth, color = :black, linestyle = :dot)

ax1.title = "γ"
ax2.title = "ωᵣ"

ax1.xlabel = "kyρ"
ax2.xlabel = "kyρ"

hlines!(ax1, [0], color = :black, linewidth=linewidth)
vlines!(ax1, [0], color = :black, linewidth=linewidth)
hlines!(ax2, [0], color = :black, linewidth=linewidth)
vlines!(ax2, [0], color = :black, linewidth=linewidth)

xlims!(ax1, 0, 0.09) 
xlims!(ax2, 0, 0.09) 

Legend(fig[1, 3],
       [gblines4, gblines5, gblines6, gblines7],
       ["β = 1.6%", "β = 1.8%", "β = 2.0%", "β = 2.2%"])

fig
