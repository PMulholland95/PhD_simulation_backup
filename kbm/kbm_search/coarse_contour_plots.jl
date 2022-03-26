using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc0 = readdlm("omtn_0_coarse_scan")
sc1 = readdlm("omt_1_omn_0_coarse_scan")
sc2 = readdlm("omt_0_omn_1_coarse_scan")
sc3 = readdlm("omt_1_omn_1_coarse_scan")

n=8

β = sc0[1:n,3];
β = convert(Array{Float64}, β);

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

ky = [0.010, 0.025, 0.050, 0.100, 0.150]


#Creating matrices
gm0 = [gam0[1:n] gam0[n+1:2n] gam0[2n+1:3n] gam0[3n+1:4n] gam0[4n+1:5n]]
gm1 = [gam1[1:n] gam1[n+1:2n] gam1[2n+1:3n] gam1[3n+1:4n] gam1[4n+1:5n]]
gm2 = [gam2[1:n] gam2[n+1:2n] gam2[2n+1:3n] gam2[3n+1:4n] gam2[4n+1:5n]]
gm3 = [gam3[1:n] gam3[n+1:2n] gam3[2n+1:3n] gam3[3n+1:4n] gam3[4n+1:5n]]

om0 = [omg0[1:n] omg0[n+1:2n] omg0[2n+1:3n] omg0[3n+1:4n] omg0[4n+1:5n]]  
om1 = [omg1[1:n] omg1[n+1:2n] omg1[2n+1:3n] omg1[3n+1:4n] omg1[4n+1:5n]]
om2 = [omg2[1:n] omg2[n+1:2n] omg2[2n+1:3n] omg2[3n+1:4n] omg2[4n+1:5n]]
om3 = [omg3[1:n] omg3[n+1:2n] omg3[2n+1:3n] omg3[3n+1:4n] omg3[4n+1:5n]]

f = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq0 = contourf!(β, ky, om0, levels = 35, colormap = :vik)

Axis(f[1, 3], title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=0 ; omnₑ=0)", xlabel = "β", ylabel = "ky")
cogam0 = contourf!(β, ky, gm0, levels = 50, colormap = :vik)

Axis(f[2, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq1 = contourf!(β, ky, om1, levels = 35, colormap = :vik)

Axis(f[2, 3], title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=1 ; omnₑ=0)", xlabel = "β", ylabel = "ky")
cogam1 = contourf!(β, ky, gm1, levels = 50, colormap = :vik)

Axis(f[3, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq2 = contourf!(β, ky, om2, levels = 35, colormap = :vik)

Axis(f[3, 3], title = "γ (omtᵢ=3 ; omnᵢ=1 ; omtₑ=0 ; omnₑ=1)", xlabel = "β", ylabel = "ky")
cogam2 = contourf!(β, ky, gm2, levels = 50, colormap = :vik)

Axis(f[4, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky") 
cofreq3 = contourf!(β, ky, om3, levels = 35, colormap = :vik)

Axis(f[4, 3], title = "γ (omtᵢ=3 ; omnᵢ=1 ; omtₑ=1 ; omnₑ=1)", xlabel = "β", ylabel = "ky")
cogam3 = contourf!(β, ky, gm3, levels = 50, colormap = :vik)


Colorbar(f[1, 2], cofreq0)

Colorbar(f[1, 4], cogam0)

Colorbar(f[2, 2], cofreq1)

Colorbar(f[2, 4], cogam1)

Colorbar(f[3, 2], cofreq2)

Colorbar(f[3, 4], cogam2)

Colorbar(f[4,2], cofreq3)

Colorbar(f[4,4], cogam3)

save("all_contour.png",f)

f
