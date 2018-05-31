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



function runModelForSobolParallel_FibrinOnly()
	allparams = readdlm("../sensitivity/sobol_paramsN1000pm50_FibrinParamsAndIC.txt", ' ', Float64);
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	outputfn="../sensitivity/sobolboundspm50percentOnlyFibrin_05_12_17.txt"
	best = mean(startingpt,1)
	@show size(allparams)
	numSets = size(allparams,1)
	paramsPerThread = numSets/nworkers()
	@sync @parallel for j in collect(1:nworkers())
	#for j in collect(1:nworkers())
		touch(string("../sensitivity/06_12_17_AUCForSobolPM50PercentN1000_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("../sensitivity/06_12_17_AUCForSobolPM50PercentN1000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			if(nworkers()==1)
				offset = 0
			else
				offset = (myid()-2)*paramsPerThread
			end
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currparams = allparams[Int((offset)+k),:]
				#need to set the non fibrin related params
				#ic at the end
				completeparams =vcat(best[1:38], currparams[1:6], best[45:47], currparams[7:36])
				#@show size(completeparams)
				fibrinIC = currparams[37:end]
				#@show size(fibrinIC)
				AUC =runModelWithParamsSetICReturnAUCFibrinIC(completeparams, fibrinIC)
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

function runModelForSobolParallel_OnlyInitialConditions_calculateCurveStats()
	#allIC = readdlm("../sensitivity/sobol_samplesOnlyIC_n2000_pm50_30_05_17.txt", ' ', Float64);
	allIC = readdlm("../sensitivity/sobol_samplesOnlyIC_n2000_pm50_02_05_18.txt", ' ', Float64);
	#startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	startingpt = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	meanparams = startingpt#mean(startingpt,1)
	params = meanparams
	tPA = 2.0
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	@show size(allIC)
	numSets = size(allIC,1)
	paramsPerThread = numSets/nworkers()
	@show paramsPerThread
	@sync @parallel for j in collect(1:nworkers())
	#for j in collect(1:nworkers())
		touch(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			if(nworkers() != 1)
				offset = (myid()-2)*paramsPerThread
			else
				offset = 0
			end
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currIC = allIC[Int((offset)+k),:]
				A,t =runModelWithParamsChangeICReturnA(params,currIC)
				AUC = calculateAUC(t,A)
				CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,t)
				write(f, string(offset+k, ",", CT, ",", CFT, ",", alpha, ",", MCF, ",", A10, ",", A20, ",", LI30, ",", LI60, ",", AUC,"\n"))
			end
		end
		close(f)
	end
end

function runModelForSobolParallel_OnlyParams_calculateCurveStats()
	#allIC = readdlm("../sensitivity/sobol_samplesOnlyIC_n2000_pm50_30_05_17.txt", ' ', Float64);
	allParams = readdlm("../sensitivity/sobol_samplesOnlyParams_n1000_pm50_02_05_18.txt", ' ', Float64);
	#startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	startingpt = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	meanparams = startingpt#mean(startingpt,1)
	params = meanparams
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	tPA = 2.0 #run this at tPA = 2 case
	@show size(allParams)
	numSets = size(allParams,1)
	paramsPerThread = numSets/nworkers()
	@show paramsPerThread
	@sync @parallel for j in collect(1:nworkers())
	#for j in collect(1:nworkers())
		touch(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyParamsN1000_", myid()-1, "_of_",nworkers(), ".txt"))
		f = open(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyParamsN1000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
		 for k in collect(1:paramsPerThread)
			if(nworkers() != 1)
				offset = (myid()-2)*paramsPerThread
			else
				offset = 0
			end
			@printf("On set %d of %d on threads %d \n", offset+k, numSets, myid())
			if(offset+k<=numSets)
				currparams = allParams[Int((offset)+k),:]
				
				A,t =runModelWithParamsReturnA(currparams,tPA)
				AUC = calculateAUC(t,A)
				CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,t)
				write(f, string(offset+k, ",", CT, ",", CFT, ",", alpha, ",", MCF, ",", A10, ",", A20, ",", LI30, ",", LI60, ",", AUC,"\n"))
			end
		end
		close(f)
	end
end

function runModelForSobolOnlyInitialConditions_calculateCurveStats_SelectIter(selIters)
	#allIC = readdlm("../sensitivity/sobol_samplesOnlyIC_n2000_pm50_30_05_17.txt", ' ', Float64);
	allIC = readdlm("../sensitivity/sobol_samplesOnlyIC_n2000_pm50_02_05_18.txt", ' ', Float64);
	#startingpt =  readdlm("../parameterEstimation/Best2PerObjectiveParameters_25_05_2017OriginalShapeFunctionOnlyFittingtPA2.txt")
	startingpt = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	meanparams = startingpt#mean(startingpt,1)
	params = meanparams
	tPA = 2.0
	#allparams = SharedArray(Float64, size(allparamslocal))
	#allparams =allparamslocal
	@show size(allIC)
	numSets = size(allIC,1)
	paramsPerThread = numSets/nworkers()
	@show paramsPerThread
	for j in selIters
	#for j in collect(1:nworkers())
		#touch(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"))
		#f = open(string("../sensitivity/02_05_18_MetricsForSobolPM50PercentOnlyICN2000_", myid()-1, "_of_",nworkers(), ".txt"), "a+")
				currIC = allIC[Int(j),:]
				A,t =runModelWithParamsChangeICReturnA(params,currIC)
				AUC = calculateAUC(t,A)
				CT,CFT,alpha,MCF,A10,A20,LI30,LI60=calculateCommonMetrics(A,t)
				@show string(j, ",", CT, ",", CFT, ",", alpha, ",", MCF, ",", A10, ",", A20, ",", LI30, ",", LI60, ",", AUC,"\n")
	end
end

