#old utilites functions
function makePredictionsFigure()
	font2 = Dict(
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>28)
	close("all")
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	ids = [3,4,9,10]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_28_11_2017.pdf")
			datasavestr = string("../generatedData/Patient", ids[j],"_tPA=", tPAs[k], "_28_11_2017.txt" )
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
			println("For patient", ids[j], "with ", tPAs[k], "tPA")
			@show CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(meanROTEM,TSIM)
			writetodisk= hcat(TSIM, transpose(meanROTEM), transpose(stdROTEM))
			writedlm(datasavestr, writetodisk)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			#plt[:subplot](size(ids,1),size(tPAs,1),counter)
			curraxis=axarr[j,k]
			plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:tick_params]("both", labelsize=18)
			curraxis[:locator_params](nbins = 4, axis = "y")
			curraxis[:locator_params](nbins = 4, axis = "x")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				curraxis[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				curraxis[:set_ylim](0, 80)
				curraxis[:set_xlim](0,TSIM[end])
			end
			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontsize=32)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
#	annotate("tPA = 0 micromolar",
#               xy=[.12;.95],
#               xycoords="figure fraction",
#               xytext=[.39,0.95],
#               textcoords="figure fraction",
#               ha="right",
#               va="top", fontsize = 24)#, family = "sans-serif")
#	annotate("tPA = 2 micromolar",
#               xy=[.12;.95],
#               xycoords="figure fraction",
#               xytext=[.85,0.95],
#               textcoords="figure fraction",
#               ha="right",
#               va="top", fontsize = 24)#, family = "sans-serif")
	axarr[1,1][:set_title]("tPA = 0 nanomolar", fontsize = 32)
	axarr[1,2][:set_title]("tPA = 2 nanomolar", fontsize =32)
	axarr[4,1][:set_xlabel]("Time, in minutes",fontsize = 32)
	axarr[4,2][:set_xlabel]("Time, in minutes",fontsize = 32)
	fig[:text](0.06, 0.5, "Clot Amplitude (mm) \n", ha="center", va="center", rotation="vertical",fontsize=40)

	figure(1)
	savefig("../figures/PredictionsFigureUsingBest2ParamSetPerObj_31_05_17OriginalShapeFunctionOnlyFittPA95Percent.pdf")
end

function makePredictionFigurePolymerizedPlatelets()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#POETs_data = "../parameterEstimation/POETS_info_28_03_18_PlateletContributionToROTEM.txt"
	POETs_data="../parameterEstimation/POETS_info_09_04_18_PlateletContributionToROTEMDiffROTEM.txt"
	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [3,4,9,10]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = 2
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_28_03_2018.pdf")
			bestparams=generateNbestPerObjective(numParamSets,ec,pc)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParams(bestparams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			curraxis=axarr[j,k]
			#plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:plot](TSIM, transpose(meanROTEM), "k")
			upper = transpose(meanROTEM+stdROTEM)
			lower = transpose(meanROTEM-stdROTEM)
			curraxis[:fill_between]((TSIM), vec(upper), vec(lower), color = ".5", alpha =.5)
			curraxis[:plot](expdata[:,1], expdata[:,2], ".k")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				curraxis[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	figure(1)
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig(string("../figures/PredictionFigureUsing",numParamSets, "ParameterSets_09_04_18_polymerizedPlateletsBest2PerObj.pdf"))
end

function makeTrainingFigureBestOveralParams()
	font2 = Dict(
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	pathToParams="../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt"
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_18_04_2017.pdf")
			datasavestr = string("../generatedData/Patient", ids[j],"_tPA=", tPAs[k], "_31_05_2017.txt" )
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicition(pathToParams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			writetodisk= hcat(TSIM, transpose(meanROTEM), transpose(stdROTEM))
			writedlm(datasavestr, writetodisk)
			@show counter
			curraxis=axarr[j,k]
			plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			if(mod(counter,2)==1)
				axis([0, TSIM[end], 0, 80])
				ylabel(string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 80])
			end
			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
#	annotate("tPA = 0 micromolar",
#               xy=[.12;.95],
#               xycoords="figure fraction",
#               xytext=[.39,0.95],
#               textcoords="figure fraction",
#               ha="right",
#               va="top", fontsize = 24, family = "sans-serif")
#	annotate("tPA = 2 micromolar",
#               xy=[.12;.95],
#               xycoords="figure fraction",
#               xytext=[.85,0.95],
#               textcoords="figure fraction",
#               ha="right",
#               va="top", fontsize = 24, family = "sans-serif")
	figure = fig
	axisarr[1,1][:set_title]("tPA = 0 nanomolar")
	axisarr[1,2][:set_title]("tPA = 2 nanomolar")

	savefig("../figures/TrainingFigureUsingBest2ParamSetPerObj_25_05_17OriginalShapeFunctionOnlyFittPA2.pdf")
end

function makeTrainingFigurePolymerizedPlatelets()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#POETs_data = "../parameterEstimation/POETS_info_28_03_18_PlateletContributionToROTEM.txt"
	POETs_data="../parameterEstimation/POETS_info_09_04_18_PlateletContributionToROTEMDiffROTEM.txt"
	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = 2
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			savestr = string("../figures/Patient", ids[j], "_tPA=", tPAs[k], "_28_03_2018.pdf")
			bestparams=generateNbestPerObjective(numParamSets,ec,pc)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParams(bestparams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			curraxis=axarr[j,k]
			#plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:plot](TSIM, transpose(meanROTEM), "k")
			upper = transpose(meanROTEM+stdROTEM)
			lower = transpose(meanROTEM-stdROTEM)
			curraxis[:fill_between]((TSIM), vec(upper), vec(lower), color = ".5", alpha =.5)
			curraxis[:plot](expdata[:,1], expdata[:,2], ".k")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				curraxis[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	figure(1)
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig(string("../figures/trainingFigureUsing",numParamSets, "ParameterSets_09_04_18_polymerizedPlateletsBest2PerObj.pdf"))
end

function makeTrainingFigure()
	font2 = Dict("family"=>"sans-serif",
	    "color"=>"black",
	    "weight"=>"normal",
	    "size"=>20)
	close("all")
	#pathToParams="parameterEstimation/Best11OverallParameters_19_04_2017.txt"
	POETs_data = "parameterEstimation/POETS_info_24_05_2017maxstep1_OriginalShapeFunction.txt"
	ec,pc,ra=parsePOETsoutput(POETs_data)
	ids = [5,6,7,8]
	tPAs = [0,2]
	close("all")
	fig,axarr = subplots(4,2,sharex="col",figsize=(15,15))
	counter = 1
	numParamSets = 15
	for j in collect(1:size(ids,1))
		for k in collect(1:size(tPAs,1))
			#savestr = string("figures/Patient", ids[j], "_tPA=", tPAs[k], "_18_04_2017.pdf")
			savestr = string("figures/Patient", ids[j], "_tPA=", tPAs[k], "_26_03_2018.pdf")
			bestparams=generateNbestGivenObjective(numParamSets,ec, pc,counter)
			alldata, meanROTEM, stdROTEM,TSIM=testROTEMPredicitionGivenParams(bestparams, ids[j], tPAs[k], savestr)
			platelets,expdata = setROTEMIC(tPAs[k], ids[j])
			@show counter
			curraxis=axarr[j,k]
			#plotAverageROTEMWDataSubplot(curraxis,TSIM,meanROTEM,stdROTEM,expdata)
			curraxis[:plot](t, transpose(meanROTEM), "k")
			upper = transpose(meanROTEM+stdROTEM)
			lower = transpose(meanROTEM-stdROTEM)
			curraxis[:fill_between]((t), vec(upper), vec(lower), color = ".5", alpha =.5)
			curraxis[:plot](expdata[:,1], expdata[:,2], ".k")
			if(mod(counter,2)==1)
				curraxis[:set_ylim](0, 70)
				curraxis[:set_xlim](0,TSIM[end])
				currasix[:set_ylabel](string("Patient ", ids[j]), fontdict=font2)
			else
				axis([0, TSIM[end], 0, 90])
			end

			if(counter==7 || counter ==8)
				xlabel("Time, in minutes", fontdict = font2)
			else
				ax =gca()
				ax[:xaxis][:set_ticklabels]([]) #remove tick labels if we're not at the bottom of a column
			end
			counter=counter+1
		end
	end
	#label columns
	annotate("tPA = 0 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.39,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")
	annotate("tPA = 2 micromolar",
               xy=[.12;.95],
               xycoords="figure fraction",
               xytext=[.85,0.95],
               textcoords="figure fraction",
               ha="right",
               va="top", fontsize = 24, family = "sans-serif")

	savefig(string("figures/trainingFigureUsing",numParamSets, "ParameterSets_22_05_17OriginalShapeFunction.pdf"))
end
