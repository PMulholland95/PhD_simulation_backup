using DelimitedFiles, CairoMakie 

sc0 = readdlm("p20sc7_scan")

n=15

β = sc0[1:n,3];
β = convert(Array{Float64}, β);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);

ky = [0.001, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


#Creating matrices
gm0 = [gam0[1:n] gam0[n+1:2n] gam0[2n+1:3n] gam0[3n+1:4n] gam0[4n+1:5n] gam0[5n+1:6n] gam0[6n+1:7n] gam0[7n+1:8n] gam0[8n+1:9n] gam0[9n+1:10n] gam0[10n+1:11n] gam0[11n+1:12n] gam0[12n+1:13n] gam0[13n+1:14n] gam0[14n+1:15n] gam0[15n+1:16n] gam0[16n+1:17n] gam0[17n+1:18n] gam0[18n+1:19n] gam0[19n+1:20n] gam0[20n+1:21n] gam0[21n+1:22n] gam0[22n+1:23n] gam0[23n+1:24n] gam0[24n+1:25n] gam0[25n+1:26n] gam0[26n+1:27n]]

om0 = [omg0[1:n] omg0[n+1:2n] omg0[2n+1:3n] omg0[3n+1:4n] omg0[4n+1:5n] omg0[5n+1:6n] omg0[6n+1:7n] omg0[7n+1:8n] omg0[8n+1:9n] omg0[9n+1:10n] omg0[10n+1:11n] omg0[11n+1:12n] omg0[12n+1:13n] omg0[13n+1:14n] omg0[14n+1:15n] omg0[15n+1:16n] omg0[16n+1:17n] omg0[17n+1:18n] omg0[18n+1:19n] omg0[19n+1:20n] omg0[20n+1:21n] omg0[21n+1:22n] omg0[22n+1:23n] omg0[23n+1:24n] omg0[24n+1:25n] omg0[25n+1:26n] omg0[26n+1:27n]]

#colormap = :gnuplot
#colormap = :jet1
#colormap = :vik
colormap = :nipy_spectral
#colormap = :plasma
#colormap = :balance
#colormap = :phase
#colormap = :matter

levels = 25
resolution = (2000, 1000)

f = Figure(resolution=resolution)

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq0 = contourf!(β, ky, om0, levels = levels, colormap = colormap)

Axis(f[1, 3], title = "γ (omtᵢ=4 ; omnᵢ=0 ; omtₑ=0 ; omnₑ=0)"
, xlabel = "β", ylabel = "ky")
cogam0 = contourf!(β, ky, gm0, levels = levels, colormap = colormap)

Colorbar(f[1, 2], cofreq0)

Colorbar(f[1, 4], cogam0)

f
