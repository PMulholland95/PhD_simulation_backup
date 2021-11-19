using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()


sc0 = readdlm("p20sc7_scan")
sc1 = readdlm("p20sc6_scan")
sca = readdlm("p24_sc0_scan")
scb = readdlm("p24_sc1_scan")

# Only use the following for p24_sc0_scan and p24_sc1_scan
#=
gama = sca[:,2];
gama = convert(Array{Float64}, gama);
omga = sca[:,2];
omga = convert(Array{Float64}, omga);  

gamb = scb[:,7];
gamb = convert(Array{Float64}, gamb);
omgb = scb[:,8];
omgb = convert(Array{Float64}, omgb);  

gam0 = [gama; gamb]
omg0 = [omga; omgb]
=#
# end of alternative approach for p24

n=15

beta = 0.002:0.002:0.03
beta = convert(Array{Float64}, beta);

#
gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);
#

ky = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

gb1a = gam0[1:n:end];
gb2a = gam0[2:n:end];
gb3a = gam0[3:n:end];
gb4a = gam0[4:n:end];
gb5a = gam0[5:n:end];
gb6a = gam0[6:n:end];
gb7a = gam0[7:n:end];
gb8a = gam0[8:n:end];
gb9a = gam0[8:n:end];
gb10a = gam0[10:n:end];
gb11a = gam0[11:n:end];
gb12a = gam0[12:n:end];
gb13a = gam0[13:n:end];
gb14a = gam0[14:n:end];
gb15a = gam0[15:n:end];

ob1a = omg0[1:n:end];
ob2a = omg0[2:n:end];
ob3a = omg0[3:n:end];
ob4a = omg0[4:n:end];
ob5a = omg0[5:n:end];
ob6a = omg0[6:n:end];
ob7a = omg0[7:n:end];
ob8a = omg0[8:n:end];
ob9a = omg0[9:n:end];
ob10a = omg0[10:n:end];
ob11a = omg0[11:n:end];
ob12a = omg0[12:n:end];
ob13a = omg0[13:n:end];
ob14a = omg0[14:n:end];
ob15a = omg0[15:n:end];


f = Figure()
fontsize_theme = Theme(fontsize = 5) 

Axis(f[1,1], title = "ωᵣ - β = 0.2%", xlabel = L"ky")
scatter!(ky, ob1a, color = :blue)

Axis(f[1,2], title = "γ - β = 0.2%", xlabel = L"ky")
scatter!(ky, gb1a, color = :red)

Axis(f[2,1], title = "ωᵣ - β = 0.4%", xlabel = L"ky")
scatter!(ky, ob2a, color = :blue)

Axis(f[2,2], title = "γ - β = 0.4%", xlabel = L"ky")
scatter!(ky, gb2a, color = :red)

Axis(f[3,1], title = "ωᵣ - β = 0.6%", xlabel = L"ky")
scatter!(ky, ob3a, color = :blue)

Axis(f[3,2], title = "γ - β = 0.6%", xlabel = L"ky")
scatter!(ky, gb3a, color = :red)

Axis(f[4,1], title = "ωᵣ - β = 0.8%", xlabel = L"ky")
scatter!(ky, ob4a, color = :blue)

Axis(f[4,2], title = "γ - β = 0.8%", xlabel = L"ky")
scatter!(ky, gb4a, color = :red)

Axis(f[1,3], title = "ωᵣ - β = 1.0%", xlabel = L"ky")
scatter!(ky, ob5a, color = :blue)

Axis(f[1,4], title = "γ - β = 1.0%", xlabel = L"ky")
scatter!(ky, gb5a, color = :red)

Axis(f[2,3], title = "ωᵣ - β = 1.2%", xlabel = L"ky")
scatter!(ky, ob6a, color = :blue)

Axis(f[2,4], title = "γ - β = 1.2%", xlabel = L"ky")
scatter!(ky, gb6a, color = :red)

Axis(f[3,3], title = "ωᵣ - β = 1.4%", xlabel = L"ky")
scatter!(ky, ob7a, color = :blue)

Axis(f[3,4], title = "γ - β = 1.4%", xlabel = L"ky")
scatter!(ky, gb7a, color = :red)

Axis(f[4,3], title = "ωᵣ - β = 1.6%", xlabel = L"ky")
scatter!(ky, ob8a, color = :blue)

Axis(f[4,4], title = "γ - β = 1.6%", xlabel = L"ky")
scatter!(ky, gb8a, color = :red)

f 

