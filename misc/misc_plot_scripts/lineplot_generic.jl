using CairoMakie
using LaTeXStrings

function plotFunction(sc::String,
		n::int,
		title::String;
		resolution=(1600,900),
		fontsize=80,
		lineswidth=10,
		)

sc0 = readdlm(sc)

beta = sc0[1:n,3]
beta = convert(Array{Float64}, beta);

kyrho = sc0[1:n:end,5];
kyrho = convert(Array{Float64}, kyrho);

gam0 = sc0[:,7];
gam0 = convert(Array{Float64}, gam0);
omg0 = sc0[:,8];
omg0 = convert(Array{Float64}, omg0);



	f = Figure()

	fontsize_theme = Theme(fontsize = fontsize)
	set_theme!(fontsize_theme)


	f
end

