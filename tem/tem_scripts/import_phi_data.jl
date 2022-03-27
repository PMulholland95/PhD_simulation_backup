# Importing phi(z) data 

using DelimitedFiles

# Extract geometry from gist: B, gyy, κ

gist = readdlm("start_data/gist_d3d_test.dat", skipstart=0);

s1 = "zprofileions_000"
s3 = ".dat"

v = []

for i in 1:9
	push!(v, string(s1, i, s3))
end

s2 = "zprofileions_00" 

for i in 10:44
	push!(v, string(s2, i, s3))
end

sc = []
p = []

for i in 1:44
	push!(sc, readdlm(string("start_data/d3d_phi_profiles/",v[i]),skipstart=7))
	push!(p, sc[i][:,3])
end

z = sc[1][:,1]

writedlm("saved_data/d3d_z.dat",z)
writedlm("saved_data/d3d_phi_z.dat",p)

