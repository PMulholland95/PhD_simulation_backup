using DelimitedFiles, CairoMakie, LaTeXStrings 

sc = readdlm("fluxspec_p42");
sc1 = readdlm("fluxspec_p02");

emk = sc[141:172,1]
emk = convert(Array{Float64}, emk); 
emq = sc[141:172,3]
emq = convert(Array{Float64}, emq);

esk = sc1[141:172,1]   
esk = convert(Array{Float64}, esk); 
esq = sc1[141:172,3]
esq = convert(Array{Float64}, esq); 

linewidth = 14
lw2 = 4
fontsize = 100
resolution = (3000,2000)

xax = 0:0.2:1.4
xax = convert(Array{Float64}, xax);
yax = 0:0.1:0.5
yax = convert(Array{Float64}, yax);

xmx1=1.24 
ymx1=0.45

fig = Figure(resolution = resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"\mathrm{k_y \rho}", ylabel = L"\mathrm{Q_i^{es} / (c_s \rho_s^2 n_{e0} T_{e0} /a^2)}", xticks = 0:0.2:1.2, xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

emql = lines!(ax1, emk, emq, color = :green, linewidth = linewidth)
esql = lines!(ax1, esk, esq, color = :blue, linewidth = linewidth) 

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth)
hlines!(ax1, [ymx1], color = :black, linewidth = linewidth)
vlines!(ax1, [xmx1], color = :black, linewidth = linewidth)
hlines!(ax1, yax, color = :black, linewidth = lw2)
vlines!(ax1, xax, color = :black, linewidth = lw2)

xlims!(ax1, 0, xmx1)
ylims!(ax1, 0, ymx1)

save("es_em_fluxspec.png", fig)

fig
