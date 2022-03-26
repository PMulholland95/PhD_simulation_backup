using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc0 = readdlm("edit_omtn_0_fine_scan")

n=11

β = sc0[1:n,3];
β = convert(Array{Float64}, β);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);

ky = [0.01, 0.02, 0.03, 0.04, 0.05]


#Creating matrices
gm0 = [gam0[1:n] gam0[n+1:2n] gam0[2n+1:3n] gam0[3n+1:4n] gam0[4n+1:5n]]

om0 = [omg0[1:n] omg0[n+1:2n] omg0[2n+1:3n] omg0[3n+1:4n] omg0[4n+1:5n]]  

f = Figure()

fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Axis(f[1, 1], title = "ωᵣ", xlabel = "β", ylabel = "ky")
cofreq0 = contourf!(β, ky, om0, levels = 100, colormap = :vik)

Axis(f[1, 3], title = "γ (omtᵢ=3 ; omnᵢ=0 ; omtₑ=0 ; omnₑ=0)"
, xlabel = "β", ylabel = "ky")
cogam0 = contourf!(β, ky, gm0, levels = 100, colormap = :vik)

Colorbar(f[1, 2], cofreq0)

Colorbar(f[1, 4], cogam0)

save("save_contour.png",f)

f
