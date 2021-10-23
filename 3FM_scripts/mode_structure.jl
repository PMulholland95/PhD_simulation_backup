using DelimitedFiles, CairoMakie, GLMakie
GLMakie.activate!()

include("convert.jl")

#mn = mode number
kxn = range(1,21, length=21)
kyn = range(1,10, length=10)

function cvint(x)
	x = convert(Array{Int}, x)
end

kxn = cvint(kxn)
kyn = cvint(kyn)

#eve = eigenvectors
eve = fm.eigenvectors


function ϕ(i,j)

	fig = Figure()

	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)

	ax1 = Axis(fig[1,1])

	mpts = fm.points[kxn[i], kyn[j]]
	mode = eve[kxn[i], kyn[j]]
	lines(ax1, mpts, map(x->abs(mode(x)), mpts))

	ax1.title = "ϕ"
	ax1.xlabel = "θ"
	
	save("single_mode.png",fig)

	fig

end



function msplots(a,b,c,d,e,f)

	fig = Figure()

	fontsize_theme = Theme(fontsize = 25)
	set_theme!(fontsize_theme)

	ax1 = Axis(fig[1,1])

	mpts1 = fm.points[kxn[a], kyn[b]]
	mode1 = eve[kxn[a], kyn[b]]

	mpts2 = fm.points[kxn[c], kyn[d]]
	mode2 = eve[kxn[c], kyn[d]]

	mpts3 = fm.points[kxn[e], kyn[f]]
	mode3 = eve[kxn[e], kyn[f]]

	ms1 = lines!(ax1, mpts1, map(x->abs(mode1(x)), mpts1), color = :red) 
	ms2 = lines!(ax1, mpts2, map(x->abs(mode2(x)), mpts2), color = :blue) 
	ms3 = lines!(ax1, mpts3, map(x->abs(mode3(x)), mpts3), color = :green) 

	ax1.title = "ϕ"
	ax1.xlabel = "θ"
	
	save("multi_mode.png",fig) 

	fig
end
