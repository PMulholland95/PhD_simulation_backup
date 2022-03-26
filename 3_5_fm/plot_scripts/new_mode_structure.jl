using CairoMakie
using LaTeXStrings

lw=12
lw2=2
xax = -15:5:15
xax = convert(Array{Float64}, xax);
yax = 0:0.5:1
yax = convert(Array{Float64}, yax);
xmx1 = 16
ymx1 = 1.1
xax2 = -25:5:25
xax2 = convert(Array{Float64}, xax2);
yax2 = 0:0.5:1
yax2 = convert(Array{Float64}, yax2);
xmx2 = 26
ymx2 = 1.1
t1 = 50
t2 = 30

function ϕ(model::PlasmaTurbulenceSaturationModel.FluidModel;
		kx=Float64,
		ky=Float64,
		resolution=(3000,2000),
		fontsize=100,
		title="Mode structure = ϕ",
		)
	
	nkx, nky = PTSM.index(kx,ky,model.spectralGrid)

	eve = model.eigenvectors
	pts = model.points
	mabs = map(x->abs(eve[nkx,nky][x]), 1:length(pts[nkx,nky])) 

	f = Figure(resolution = resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	ax1 =   Axis(f[1,1], 
		xlabel = L"\theta / \pi",
		ylabel = L"\phi", 
		xminorticksvisible = true, 
		yminorticksvisible = true, 
		xticksize = t1, 
		yticksize = t1, 
		xminorticksize = t2, 
		yminorticksize = t2)

	lines!(pts[nkx,nky], mabs, color = :blue, linewidth = lw)

	hlines!(ax1, [0], color = :black, linewidth = lw)
	hlines!(ax1, [ymx1], color = :black, linewidth = lw)
	vlines!(ax1, [xmx1], color = :black, linewidth = lw)
	vlines!(ax1, [-xmx1], color = :black, linewidth = lw)
	hlines!(ax1, yax, color = :black, linewidth = lw2)
	vlines!(ax1, xax, color = :black, linewidth = lw2)

	xlims!(ax1, -xmx1, xmx1)
	ylims!(ax1, 0, ymx1)

	save("single_mode_struct.png", f)

	f

end



function msplots(model::PlasmaTurbulenceSaturationModel.FluidModel;
		kxa=Float64,
		kya=Float64,
		kxb=Float64,
		kyb=Float64,
		kxc=Float64,
		kyc=Float64,
		resolution=(3000,2000),
		fontsize=100,
		title="Mode structure = ϕ",
		)

	nkxa, nkya = PTSM.index(kxa,kya,model.spectralGrid)
	nkxb, nkyb = PTSM.index(kxb,kyb,model.spectralGrid)
	nkxc, nkyc = PTSM.index(kxc,kyc,model.spectralGrid)

	eve = model.eigenvectors
	pts = model.points
	#evec = circshift(eve, (10,0))
	#ptsc = circshift(pts, (10,0))

	mabsa = map(x->abs(eve[nkxa,nkya][x]), 1:length(pts[nkxa,nkya])) 
	mabsb = map(x->abs(eve[nkxb,nkyb][x]), 1:length(pts[nkxb,nkyb])) 
	mabsc = map(x->abs(eve[nkxc,nkyc][x]), 1:length(pts[nkxc,nkyc]))

	f = Figure(resolution=resolution)

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)


	ax1 =   Axis(f[1,1], 
		xlabel = L"\theta / \pi",
		ylabel = L"\phi", 
		xminorticksvisible = true, 
		yminorticksvisible = true, 
		xticksize = t1, 
		yticksize = t1, 
		xminorticksize = t2, 
		yminorticksize = t2)

	lines!(pts[nkxa,nkya], mabsa, color = :red, linewidth=lw)
	lines!(pts[nkxb,nkyb], mabsb, color = :blue, linewidth=lw)
	lines!(pts[nkxc,nkyc], mabsc, color = :green, linewidth=lw)
	
	hlines!(ax1, [0], color = :black, linewidth = lw)
	hlines!(ax1, [ymx2], color = :black, linewidth = lw)
	vlines!(ax1, [xmx2], color = :black, linewidth = lw)
	vlines!(ax1, [-xmx2], color = :black, linewidth = lw)
	hlines!(ax1, yax2, color = :black, linewidth = lw2)
	vlines!(ax1, xax2, color = :black, linewidth = lw2)

	xlims!(ax1, -xmx2, xmx2)
	ylims!(ax1, 0, ymx1)

	save("multiple_mode_struct.png", f)

	f
end
