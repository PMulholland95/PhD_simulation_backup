using DelimitedFiles, CairoMakie 

sc0 = readdlm("full_consist_scan")

n=35
resolution = (3000,2000) 

β = sc0[1:n,3];
β = convert(Array{Float64}, β);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);

ky = 0.05:0.05:0.8


#Creating matrices
gm0 = [gam0[1:n] gam0[n+1:2n] gam0[2n+1:3n] gam0[3n+1:4n] gam0[4n+1:5n] gam0[5n+1:6n] gam0[6n+1:7n] gam0[7n+1:8n] gam0[8n+1:9n] gam0[9n+1:10n] gam0[10n+1:11n] gam0[11n+1:12n] gam0[12n+1:13n] gam0[13n+1:14n] gam0[14n+1:15n] gam0[15n+1:16n]]

om0 = [omg0[1:n] omg0[n+1:2n] omg0[2n+1:3n] omg0[3n+1:4n] omg0[4n+1:5n] omg0[5n+1:6n] omg0[6n+1:7n] omg0[7n+1:8n] omg0[8n+1:9n] omg0[9n+1:10n] omg0[10n+1:11n] omg0[11n+1:12n] omg0[12n+1:13n] omg0[13n+1:14n] omg0[14n+1:15n] omg0[15n+1:16n]]

f = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = 80)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq0 = contourf!(β, ky, om0, levels = 50, colormap = :nipy_spectral)

Axis(f[1, 3], title = "γ", xlabel = "β", ylabel = "ky")
cogam0 = contourf!(β, ky, gm0, levels = 50, colormap = :nipy_spectral)

Colorbar(f[1, 2], cofreq0)

Colorbar(f[1, 4], cogam0)

save("consistent_contour.png",f)

f
