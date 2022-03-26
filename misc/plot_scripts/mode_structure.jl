using CairoMakie
using LaTeXStrings


function ϕ(model::PlasmaTurbulenceSaturationModel.FluidModel;
		kx=Float64,
		ky=Float64,
		resolution=(1600,900),
		fontsize=25,
		title="Mode structure = ϕ",
		)
	
	nkx, nky = PTSM.index(kx,ky,model.spectralGrid)

	eve = model.eigenvectors
	pts = model.points
	mabs = map(x->abs(eve[nkx,nky][x]), 1:length(pts[nkx,nky])) 

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1,1], title = title, xlabel = "θ")
	lines!(pts[nkx,nky], mabs)

	f

end



function msplots(model::PlasmaTurbulenceSaturationModel.FluidModel;
		kxa=Float64,
		kya=Float64,
		kxb=Float64,
		kyb=Float64,
		kxc=Float64,
		kyc=Float64,
		resolution=(1600,900),
		fontsize=25,
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

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1,1], title = title, xlabel = "θ")
	lines!(pts[nkxa,nkya], mabsa, color = :red)
	lines!(pts[nkxb,nkyb], mabsb, color = :blue)
	lines!(pts[nkxc,nkyc], mabsc, color = :green)
	
	f
end
