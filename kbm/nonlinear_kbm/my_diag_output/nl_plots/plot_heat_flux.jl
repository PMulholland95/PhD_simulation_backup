using DelimitedFiles, CairoMakie, LaTeXStrings 

sc = readdlm("nrg_p42"; skipstart=6);
sc1 = readdlm("nrg_p02"; skipstart=6);

emt = sc[:,1]
emq = sc[:,3]

est = sc1[:,1]   
esq = sc1[:,3]

linewidth = 14
lw2 = 4
fontsize = 100
resolution = (3000,2000)

xax = 0:250:length(est)
xax = convert(Array{Float64}, xax);
yax = 0:0.5:4
yax = convert(Array{Float64}, yax);

xmx1=1450 
ymx1=4.2

fig = Figure(resolution = resolution)

fontsize_theme = Theme(fontsize = fontsize)
set_theme!(fontsize_theme)

ax1 = Axis(fig[1,1], xlabel = L"\mathrm{t / (L_{ref}/c_s)}", ylabel = L"\mathrm{Q_i^{es} / (c_s \rho_s^2 n_{e0} T_{e0} /a^2)}", xticks = 0:250:1500, xminorticksvisible = true, yminorticksvisible = true, xticksize = 40, yticksize = 40, xminorticksize = 20, yminorticksize = 20)

emql = lines!(ax1, emt, emq, color = :green, linewidth = linewidth)
esql = lines!(ax1, est, esq, color = :blue, linewidth = linewidth) 

hlines!(ax1, [0], color = :black, linewidth = linewidth) 
vlines!(ax1, [0], color = :black, linewidth = linewidth)
hlines!(ax1, [ymx1], color = :black, linewidth = linewidth)
vlines!(ax1, [xmx1], color = :black, linewidth = linewidth)
hlines!(ax1, yax, color = :black, linewidth = lw2)
vlines!(ax1, xax, color = :black, linewidth = lw2)

xlims!(ax1, 0, xmx1)
ylims!(ax1, 0, ymx1)

save("es_em_heat_flux.png", fig)

fig
