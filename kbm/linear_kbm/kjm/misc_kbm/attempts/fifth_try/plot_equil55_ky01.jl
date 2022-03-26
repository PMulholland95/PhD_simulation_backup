using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc = readdlm("beta_equil55_ky01");

beta = sc[:,3];
beta = convert(Array{Float64}, beta);

gam = sc[:,5];
gam = convert(Array{Float64}, gam);

omg = sc[:,6];
omg = convert(Array{Float64}, omg);

# We now want to gather together gamma values for range of kyrho, and fixed beta
# In the end, only want to plot: gamma|max vs. beta 

#Makes array (gambet0) of every 6th value in original array (gam) starting from first value (beta=0)
# This gathers all gammas with beta=0 for all kyrho values in scan

# plot(kyrho, gambet0)
fig = Figure()

ax1 = Axis(fig[1,1])
ax2 = Axis(fig[1,2])

lines!(ax1, beta, gam, color = :red)
lines!(ax2, beta, omg, color = :blue)

ax1.title = "γ [kyρ = 0.1]"
ax2.title = "ω [kyρ = 0.1]"

ax1.xlabel = "β"
ax2.xlabel = "β"

xlims!(ax1, 0.03, 0.036)
xlims!(ax2, 0.03, 0.036)
ylims!(ax1, 0, 0.1)
ylims!(ax2, 0, 0.1)

save("beta_equil55_ky01.png",fig)

fig
