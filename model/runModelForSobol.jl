include("runModel.jl")

function runModelForSobol()
	allparams = readdlm("sensitivity/paramsplusorminus50percentN5000.txt", ' ', Float64);
	@show size(allparams)
	numSets = size(allparams,1)
	f = open("sensitivity/AUCForSobolPM50PercentN1000.txt", "a+")
	for j in collect(1:numSets)
		@printf("On set %d of %d\n", j, numSets)
		currparams = allparams[j,:]
		AUC =runModelWithParamsReturnAUC(currparams)
		write(f, string(AUC, "\n"))
	end
	close(f)
end

function runModelForSobolParallel()
	allparams = readdlm("/..sensitivity/sobol_samplesOnlyParams_n500_pm50_30_05_17.txt", ' ', Float64);
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	@show size(allparams)
	numSets = size(allparams,1)
	paramsPerThread = numSets/nworkers()
	@sync @parallel for j in collect(1:nworkers())
		touch(string("../sensitivity/AUCForSobolPM50PercentN500OnlyParams_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("../sensitivity/AUCForSobolPM50PercentN500OnlyParams_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			offset = (myid()-2)*paramsPerThread
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currparams = allparams[Int((offset)+k),:]
				AUC =runModelWithParamsReturnAUC(currparams,2)
				#@show AUC
				write(f, string(offset+k, ",", AUC, "\n"))
			end
		end
		close(f)
	end
end

function runModelForSobolParallel_IncludingInitialConditions()
	allparams = readdlm("../sensitivity/sobol_paramsN5000pm50_ICandParams_05_31_17.txt", ' ', Float64);
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	@show size(allparams)
	numSets = size(allparams,1)
	paramsPerThread = numSets/nworkers()
	@sync @parallel for j in collect(1:nworkers())
		touch(string("../sensitivity/05_31_17_AUCForSobolPM50PercentN5000_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("../sensitivity/05_31_17_AUCForSobolPM50PercentN5000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			offset = (myid()-2)*paramsPerThread
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currparams = allparams[Int((offset)+k),:]
				AUC =runModelWithParamsSetICReturnAUC(currparams)
				#@show AUC
				write(f, string(offset+k, ",", AUC, "\n"))
			end
		end
		close(f)
	end
end

function runModelForSobolParallel_OnlyInitialConditions()
	allIC = readdlm("sensitivity/sobol_samplesOnlyIC_n2000_pm50_30_05_17.txt", ' ', Float64);
	startingpt =  readdlm("parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	meanparams = mean(startingpt,1)
	params = meanparams
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	@show size(allIC)
	numSets = size(allIC,1)
	paramsPerThread = numSets/nworkers()
	@sync @parallel for j in collect(1:nworkers())
		touch(string("sensitivity/05_30_17_AUCForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("sensitivity/05_30_17_AUCForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			offset = (myid()-2)*paramsPerThread
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currIC = allIC[Int((offset)+k),:]
				AUC =runModelWithParamsChangeICReturnAUC(params,currIC)
				#@show AUC
				write(f, string(offset+k, ",", AUC, "\n"))
			end
		end
		close(f)
	end
end
