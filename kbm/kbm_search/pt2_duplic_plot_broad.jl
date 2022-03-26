using DelimitedFiles, CairoMakie, GLMakie 
GLMakie.activate!()

sc0 = readdlm("p24_sc0_scan")

n=15

beta = 0.002:0.002:0.03
beta = convert(Array{Float64}, beta);

gam0 = sc0[:,2];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,3];
omg0 = convert(Array{Float64}, omg0);

ky = [0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.4, 0.5]

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
fontsize_theme = Theme(fontsize = 30) 

Axis(f[1,1], title = "ωᵣ - β = 1.8%", xlabel = L"ky")
scatter!(ky, ob9a, color = :blue)

Axis(f[1,2], title = "γ - β = 1.8%", xlabel = L"ky")
scatter!(ky, gb9a, color = :red)

Axis(f[2,1], title = "ωᵣ - β = 2.0%", xlabel = L"ky")
scatter!(ky, ob10a, color = :blue)

Axis(f[2,2], title = "γ - β = 2.0%", xlabel = L"ky")
scatter!(ky, gb10a, color = :red)

Axis(f[3,1], title = "ωᵣ - β = 2.2%", xlabel = L"ky")
scatter!(ky, ob11a, color = :blue)

Axis(f[3,2], title = "γ - β = 2.2%", xlabel = L"ky")
scatter!(ky, gb11a, color = :red)

Axis(f[4,1], title = "ωᵣ - β = 2.4%", xlabel = L"ky")
scatter!(ky, ob12a, color = :blue)

Axis(f[4,2], title = "γ - β = 2.4%", xlabel = L"ky")
scatter!(ky, gb12a, color = :red)

Axis(f[1,3], title = "ωᵣ - β = 2.6%", xlabel = L"ky")
scatter!(ky, ob13a, color = :blue)

Axis(f[1,4], title = "γ - β = 2.6%", xlabel = L"ky")
scatter!(ky, gb13a, color = :red)

Axis(f[2,3], title = "ωᵣ - β = 2.8%", xlabel = L"ky")
scatter!(ky, ob14a, color = :blue)

Axis(f[2,4], title = "γ - β = 2.8%", xlabel = L"ky")
scatter!(ky, gb14a, color = :red)

Axis(f[3,3], title = "ωᵣ - β = 3.0%", xlabel = L"ky")
scatter!(ky, ob15a, color = :blue)

Axis(f[3,4], title = "γ - β = 3.0%", xlabel = L"ky")
scatter!(ky, gb15a, color = :red)

f
