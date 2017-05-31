using PyPlot
using PyCall
PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Helvetica"]

function makeSobolGraphNoAnnotations()
	font1 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>24)

	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>12)
	close("all")
	numparams = 77+22
	fig=figure(figsize=[35,15])
	plt[:tight_layout]=true
	data = readdlm("../sensitivity/soboloutputpm50percentICandParams_05_31_17_N1000.txt")
	topHalf = data[1:numparams+1, :]
	@show topHalf
	usefulData = topHalf[2:end, :]
	positions = collect(0:numparams-1)
	bar(positions, usefulData[:,4],color = "k", yerr=usefulData[:,5], align="center")
	ax = gca()
	ax[:xaxis][:set_ticks](positions)
	ylabel("Total Order Sensitivity Indicies", fontdict=font1)
	axis("tight")
	axis([0,numparams,0,.4])
	ax[:tick_params](labelsize=20)
	#lines and label for kinetic parameters
	ax[:xaxis][:set_ticklabels](usefulData[:,1], rotation = 90, fontsize = 14)
	savefig("../sensitivity/SobolTotalOrderN1000_05_31_17.pdf", bbox_inches = "tight")
end

function makeSobolGraphNoAnnotationsOnlyIC()
	font1 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>24)

	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>12)
	close("all")
	numparams =22
	fig=figure(figsize=[25,15])
	data = readdlm("sensitivity/soboloutputpm50percentOnlyIC_05_30_17N2000.txt")
	topHalf = data[1:numparams+1, :]
	@show topHalf
	usefulData = topHalf[2:end, :]
	positions = collect(1:numparams)
	bar(positions, usefulData[:,4],color = "k", yerr=usefulData[:,5], align="center")
	ax = gca()
	ax[:xaxis][:set_ticks](positions)
	ylabel("Total Order Sensitivity Indicies", fontdict=font1)
	axis("tight")
	axis([0,numparams,0,1])
	ax[:tick_params](labelsize=20)
	#lines and label for kinetic parameters
	ax[:xaxis][:set_ticklabels](usefulData[:,1], rotation = 50, fontsize = 8)
	savefig("sensitivity/SobolTotalOrderOnlyICN52000_05_30_17.pdf")
end

function makeSobolGraph()
	font1 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>24)

	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>12)
	close("all")
	numparams = 77
	fig=figure(figsize=[25,15])
	data = readdlm("sensitivity/SobolOputput_pm50_N5000_04_24_2017.txt")
	topHalf = data[1:numparams+1, :]
	@show topHalf
	usefulData = topHalf[2:end, :]
	positions = collect(0:numparams-1)
	bar(positions, usefulData[:,4],color = "k", yerr=usefulData[:,5], align="center")
	ax = gca()
	ax[:xaxis][:set_ticks](positions)
	ylabel("Total Order Sensitivity Indicies", fontdict=font1)
	axis("tight")
	axis([0,numparams,0,1])
	ax[:tick_params](labelsize=20)
	#lines and label for kinetic parameters
	ax[:xaxis][:set_ticklabels]([], rotation = 80, fontsize = 5)
	#ax[:xaxis][:set_ticklabels](usefulData[:,1], rotation = 80, fontsize = 5)
	annotate("",
		xy=[.125;.83],
		xycoords="figure fraction",
		xytext=[.3,0.83],
		textcoords="figure fraction",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))
	annotate("Thombin Generation Kinetic Parameters",
		xy=[.12;.81],
		xycoords="figure fraction",
		xytext=[.29,0.81],
		textcoords="figure fraction",
		ha="right",
		va="top")

	annotate("",
		xy=[.3;.72],
		xycoords="figure fraction",
		xytext=[.5,0.72],
		textcoords="figure fraction",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Thrombin Generation Control Parameters",
		xy=[.3;.71],
		xycoords="figure fraction",
		xytext=[.475,0.71],
		textcoords="figure fraction",
		ha="right",
		va="top")

	annotate("",
		xy=[.505;.83],
		xycoords="figure fraction",
		xytext=[.59,0.83],
		textcoords="figure fraction",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Platelet Parameters",
		xy=[.3;.82],
		xycoords="figure fraction",
		xytext=[.575,0.82],
		textcoords="figure fraction",
		ha="right",
		va="top")

	annotate("",
		xy=[.595;.64],
		xycoords="figure fraction",
		xytext=[.81,0.64],
		textcoords="figure fraction",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Fibrin Kinetic Parameters",
		xy=[.65;.63],
		xycoords="figure fraction",
		xytext=[.75,0.63],
		textcoords="figure fraction",
		ha="right",
		va="top")

		annotate("",
		xy=[.815;.85],
		xycoords="figure fraction",
		xytext=[.9,0.85],
		textcoords="figure fraction",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Fibrin Control Parameters",
		xy=[.85;.82],
		xycoords="figure fraction",
		xytext=[.9,0.82],
		textcoords="figure fraction",
		ha="right",
		va="top")
	#label columns of interest
	annotate("k_inhibition_ATIII",
		xy=[8;.35],# Arrow tip
		xycoords="data", # Coordinates in in "data" units
		xytext=[9;.45], # Text offset from tip
		textcoords="data",
		ha="right",
		va="top",
		arrowprops=Dict("facecolor"=>"black", "width"=>.5, "headwidth"=>3)) #

	annotate("Trigger Control Parameters",
		xy=[18.5,0.05],
		xycoords="data",
		xytext=[20.5,0.15],
		textcoords="data",
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "width"=>.5))

	annotate("k_cat_fibrinogen",
		xy=[47;.33],# Arrow tip
		xycoords="data", # Coordinates in in "data" units
		xytext=[46;.45], # Text offset from tip
		textcoords="data",
		ha="right",
		va="top",
		arrowprops=Dict("facecolor"=>"black", "width"=>.5, "headwidth"=>3))

	savefig("sensitivity/SobolTotalOrderN5000.pdf")
	
end
