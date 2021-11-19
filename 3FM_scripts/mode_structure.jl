using CairoMakie
using LaTeXStrings


function ϕ(model::PlasmaTurbulenceSaturationModel.FluidModel;
		nkx=Int,
		nky=Int,
		resolution=(1600,900),
		fontsize=25,
		title="Mode structure = ϕ",
		)
	
	eve = model.eigenvectors
	pts = model.points
	evec = circshift(eve, (10,0))
	ptsc = circshift(pts, (10,0))
	mabs = map(x->abs(evec[nkx,nky](x)), ptsc[nkx,nky]) 

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1,1], title = title, xlabel = "θ")
	lines!(ptsc[nkx,nky], mabs)

	f

end



function msplots(model::PlasmaTurbulenceSaturationModel.FluidModel;
		nkxa=Int,
		nkya=Int,
		nkxb=Int,
		nkyb=Int,
		nkxc=Int,
		nkyc=Int,
		resolution=(1600,900),
		fontsize=25,
		title="Mode structure = ϕ",
		)

	
	eve = model.eigenvectors
	pts = model.points
	evec = circshift(eve, (10,0))
	ptsc = circshift(pts, (10,0))

	mabsa = map(x->abs(evec[nkxa,nkya](x)), ptsc[nkxa,nkya]) 
	mabsb = map(x->abs(evec[nkxb,nkyb](x)), ptsc[nkxb,nkyb]) 
	mabsc = map(x->abs(evec[nkxc,nkyc](x)), ptsc[nkxc,nkyc]) 

	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)

	Axis(f[1,1], title = title, xlabel = "θ")
	lines!(ptsc[nkxa,nkya], mabsa, color = :red)
	lines!(ptsc[nkxb,nkyb], mabsb, color = :blue)
	lines!(ptsc[nkxc,nkyc], mabsc, color = :green)
	
	f
end
