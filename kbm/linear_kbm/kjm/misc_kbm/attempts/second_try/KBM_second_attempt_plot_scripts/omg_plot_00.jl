using DelimitedFiles, Makie, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("00scan");

beta = sc[1:6,3];
beta = convert(Array{Float64}, beta);

#gam is actually omega here (mode frequency)
gam = sc[:,8];
gam = convert(Array{Float64}, gam);

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

# plot(kyrho, gambet0)
CairoMakie.lines(kyrho, gambet0, color = :red, label = "beta=0")
CairoMakie.lines!(kyrho, gambet1, color = :orange, label = "beta=0.001")
CairoMakie.lines!(kyrho, gambet2, color = :yellow, label = "beta=0.002")
CairoMakie.lines!(kyrho, gambet3, color = :green, label = "beta=0.003")
CairoMakie.lines!(kyrho, gambet4, color = :blue, label = "beta=0.004")
CairoMakie.lines!(kyrho, gambet5, color = :purple, label = "beta=0.005")
CairoMakie.current_figure()


