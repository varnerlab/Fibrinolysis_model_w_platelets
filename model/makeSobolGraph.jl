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
	data = readdlm("../sensitivity/SobolOutpt_05_31_17_AUCForSobolParamsAndICPM50PercentN5000.txt")
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
	savefig("../sensitivity/SobolTotalOrderN5000_05_31_17.pdf", bbox_inches = "tight")
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
	font1 = Dict(
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>40)

	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>12)
	close("all")
	numparams = 77
	fig=figure(figsize=[25,15])
	data = readdlm("../sensitivity/SobolOputput_pm50_N5000_04_24_2017.txt")
	topHalf = data[1:numparams+1, :]
	@show topHalf
	usefulData = topHalf[2:end, :]
	positions = collect(0:numparams-1)
	bar(positions, usefulData[:,4],color = "k", yerr=usefulData[:,5], align="center", ecolor = "grey")
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
		xy=[0;-.02],
		xycoords=("data","axes fraction"),
		xytext=[17.5;-0.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))
	annotate("Thombin Kinetic    ",
		fontsize=30,
		xy=[0;-.03],
		xycoords=("data","axes fraction"),
		xytext=[17;-.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")

	annotate("",
		xy=[17.5;-.02],
		xycoords=("data","axes fraction"),
		xytext=[38,-.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Thrombin Control    ",
		fontsize=30,
		xy=[17.5,-.03],
		xycoords=("data","axes fraction"),
		xytext=[38,-0.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")

	annotate("",
		xy=[38.0;-.02],
		xycoords=("data","axes fraction"),
		xytext=[44,-.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Platelet  ",
		fontsize=24,
		xy=[38.25;-0.03],
		xycoords=("data","axes fraction"),
		xytext=[44,-.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")

	annotate("",
		xy=[44;-.02],
		xycoords=("data","axes fraction"),
		xytext=[47,-.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Timing",
		fontsize=20,
		xy=[44;-0.03],
		xycoords=("data","axes fraction"),
		xytext=[47,-.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")
	annotate("",
		xy=[47;-.02],
		xycoords=("data","axes fraction"),
		xytext=[69,-.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Fibrin Kinetic             ",
		fontsize=30,
		xy=[47;-.03],
		xycoords=("data","axes fraction"),
		xytext=[69,-.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")

		annotate("",
		xy=[69;-.02],
		xycoords=("data","axes fraction"),
		xytext=[77,-.02],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top",
		arrowprops = Dict("facecolor"=> "black", "arrowstyle" => "-" ))

	annotate("Fibrin  \nControl ",
		fontsize =30,
		xy=[69;-.03],
		xycoords=("data","axes fraction"),
		xytext=[77,-.03],
		textcoords=("data","axes fraction"),
		ha="right",
		va="top")
	#label columns of interest
	annotate("k_inhibition_ATIII",
		xy=[8;.35],# Arrow tip
		xycoords="data", # Coordinates in in "data" units
		xytext=[9;.45], # Text offset from tip
		textcoords="data",
		ha="center",
		va="top",
		fontsize=30,
		arrowprops=Dict("facecolor"=>"black", "width"=>.5, "headwidth"=>3)) #

	annotate("Trigger Control Parameters",
		xy=[18.5,0.05],
		xycoords="data",
		xytext=[20.5,0.15],
		textcoords="data",
		ha="center",
		va="top",
		fontsize=30,
		arrowprops = Dict("facecolor"=> "black", "width"=>.5))

	annotate("k_cat_fibrinogen",
		xy=[47;.33],# Arrow tip
		xycoords="data", # Coordinates in in "data" units
		xytext=[46;.45], # Text offset from tip
		textcoords="data",
		ha="center",
		va="top",
		fontsize=30,
		arrowprops=Dict("facecolor"=>"black", "width"=>.5, "headwidth"=>3))
	xlabel("\nParameter Index", fontsize = 40)

	savefig("../sensitivity/SobolTotalOrderN5000_24_08_17.pdf")
	
end
