using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc0 = readdlm("p22sc3_scan")

gam0 = sc0[:,5];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,6];
omg0 = convert(Array{Float64}, omg0);

ky = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


f = Figure()
fontsize_theme = Theme(fontsize = 100) 

Axis(f[1,1], xlabel = L"ky", ylabel = "ωᵣ")
scatter!(ky, omg0, color = :blue, markersize = 25, marker = '△')

Axis(f[1,2], xlabel = L"ky", ylabel = "γ")
scatter!(ky, gam0, color = :red, markersize = 25, marker = '△')

f
