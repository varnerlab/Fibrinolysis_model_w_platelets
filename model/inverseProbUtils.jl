#utils for analyzing solution from solving inverse problem
using PyPlot
using DelimitedFiles
using Statistics
using Random
using Distributions

include("runModel.jl")

function sampleSpace(lower,upper,cond)
	res = zeros(1)
	if(size(cond,1)==1)
		return rand(Normal(cond, (upper*cond-lower*cond)/6),1)[1]
	else
		for j = 1:maximum(size(cond))
			#@show (upper*cond[j]-lower*cond[j])
			if(cond[j]==0)
				temp = rand( Truncated(Normal(0, .005), 0, .01))
			else
				temp = rand(Normal(cond[j], (upper*cond[j]-lower*cond[j])/8),1)[1]
			end
			#@show temp
			push!(res, temp)
		end
	end
	#@show cond, res
	return res[2:end] #chop off first element
end

function testIfPhysical(params,genIC,genExp,genPlatelets)
	T,R=runModelWithParamsChangeICReturnA(params,genIC,genExp,genPlatelets)
	#plot(T,R)
	#check to make sure solver finished
	if(T[end]<59)
		println("Didn't finish solving")
		isPhysical=false
		return isPhysical
	end

	#describe the curve we've created'
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	target_LI30 = calculateLysisAtTime(R,T,30.0)
	target_AUC = calculateAUC(R,T)


	temp_target=[target_CT, target_CFT, target_alpha, target_MCF,target_MaximumLysis, target_AUC, target_LI30]
	@show temp_target
	if(true in (temp_target .<0))
		isPhysical = false
	else
		isPhysical = true
	end
	#based on ""Normal range values for thromboelastography in healthy adult volunteers""
	#want MCF between 50 and 70
	#MA = MCF
	#R = CT
	#K = CFT
	#if(target_MCF<49.7 || target_MCF>72.7)
	if(target_MCF<40.0 || target_MCF>72.7) #looser contstraints
		isPhysical=false
	end
	if(target_CT<3.8*60  || target_CT>9.8*60)
		isPhysical=false
	end
	if(target_CFT<.7*60  || target_CFT>3.4*60)
		isPhysical=false
	end
#	if(target_alpha<47.8  || target_CFT>77.7)
#		isPhysical=false
#	end

	return isPhysical, temp_target
end

function objective_f(params::Vector, grad::Vector)
	curr_exp = params[1:7]
	curr_ICs = params[8:end-1]
	curr_mx = params[end]
	tPA = curr_ICs[12]
	dilution_factor = 1.0
	T,X,R = solve_model_with_parameters_and_IC(parameter_array, dilution_factor, tPA, curr_mx,curr_ICs, curr_exp)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	#if any of the metrics are negative, something went wrong, penalize our parameters
	if(true in (metrics .<0)) #if any of our metrics are negative
		calc_obj = 1E8
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:maximum(size(metrics))
			sum = sum+(metrics[j]-target[j])^2
		end
		@show params
		@show metrics
		calc_obj = sqrt(sum)
	end
	@show calc_obj
	return calc_obj

end

function objective_five_metrics(params::Vector, grad::Vector)
	curr_exp = params[1:7]
	curr_ICs = params[8:end-1]
	curr_mx = params[end]
	tPA = curr_ICs[12]
	dilution_factor =  0.79
	T,X,R = solve_model_with_parameters_and_IC(parameter_array, dilution_factor, tPA, curr_mx,curr_ICs, curr_exp)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	maximumLysis = calculateLysisAtTime(R,T,120.0)
	sel_metrics = [CT,CFT,alpha,MCF,maximumLysis]
	#@show abs(sel_metrics-sel_target)
	#if any of the metrics are negative, something went wrong, penalize our parameters
	if(true in (metrics .<0)) #if any of our metrics are negative
		calc_obj = 1E8
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:maximum(size(sel_metrics))
			sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2
		end
		#@show params
		#@show metrics
		calc_obj = sqrt(sum)
	end
	#@show calc_obj
	return calc_obj

end

function objective_five_metrics_weighted(params::Vector, grad::Vector)
	curr_exp = params[1:7]
	curr_ICs = params[8:end-1]
	curr_mx = params[end]
	tPA = curr_ICs[12]
	dilution_factor =  0.79
	T,X,R = solve_model_with_parameters_and_IC(parameter_array, dilution_factor, tPA, curr_mx,curr_ICs, curr_exp)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	maximumLysis = calculateLysisAtTime(R,T,120.0)
	sel_metrics = [CT,CFT,alpha,MCF,maximumLysis]
	#@show abs(sel_metrics-sel_target)
	#if any of the metrics are negative, something went wrong, penalize our parameters
	if(true in (metrics .<0)) #if any of our metrics are negative
		calc_obj = 1E8
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:maximum(size(sel_metrics))
			sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
		end
		#@show params
		#@show metrics
		calc_obj = sqrt(sum)
	end
	#@show calc_obj
	return calc_obj

end

function objective_six_metrics_weighted(params::Vector, grad::Vector)
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	T,R =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	maximumLysis = calculateLysisAtTime(R,T,120.0)
	AUC = calculateAUC(R,T)
	sel_metrics = [CT,CFT,alpha,MCF,maximumLysis,AUC]
	#@show abs(sel_metrics-sel_target)
	#if any of the metrics are negative, something went wrong, penalize our parameters
	calc_obj = 0.0
	if(true in (sel_metrics .<0)) #if any of our metrics are negative
		#count how many are negative-the more negative, the worse
		for j in 1:maximum(size(sel_metrics))
			if(sel_metrics[j]<0)
				calc_obj=calc_obj+1E4
			end
		end
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:maximum(size(sel_metrics))
			sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
		end
		calc_obj = sqrt(sum)
	end
	@show sel_metrics, calc_obj
	return calc_obj

end

function objective_six_metrics_weighted(params)
	#@show params
	#check for negativity in params. If so, heavily penalize
	if(true in (params .<0))
		calc_obj=1E8
		return calc_obj
	end
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	@time T,R =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	#experiment run for 60 mins =1 hours
	maximumLysis = calculateLysisAtTime(R,T,60.0)
	AUC = calculateAUC(R,T)
	sel_metrics = [CT,CFT,alpha,MCF,maximumLysis,AUC]
	#@show abs(sel_metrics-sel_target)
	#if any of the metrics are negative, something went wrong, penalize our parameters
	calc_obj = 0.0
	if(true in (sel_metrics .<0)) #if any of our metrics are negative
		#count how many are negative-the more negative, the worse
		for j in 1:maximum(size(sel_metrics))
			if(sel_metrics[j]<0)
				calc_obj=calc_obj+1E4
			end
		end
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:maximum(size(sel_metrics))
			sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
		end
		calc_obj = sqrt(sum)
	end
	@show sel_metrics, calc_obj
	return calc_obj

end

function plotAllCurveSameProb()
	close("all")
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	kin_params=mean(allp, dims=1)
	originalIC = readdlm("../solveInverseProb/Master_ics_to_match_13_03_19_.txt", ',')
	curr_exp = originalIC[1:8]
	curr_ICs = originalIC[9:end-1]
	curr_platelets = originalIC[end]
	tPA = originalIC[12]
	T,R =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	figure(figsize=[15,15])
	plot(T,R, "k", linewidth=2)
	
	numSims = 1
	solved=[collect(1:6);collect(201:206)]
	allCurves = zeros(size(solved,1), size(R,1))
	count = 1
	for k in solved
		currEstICs = readdlm(string("../solveInverseProb/foundIcs_13_03_19_", k, "solvingSameProb.txt"))
		curr_exp = currEstICs[1:8]
		curr_ICs = currEstICs[9:end-1]
		curr_platelets = currEstICs[end]
		tPA = currEstICs[12]
		T,R=runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
		plot(T,R, color = "gray", alpha =.5, linewidth=.2)
		@show size(allCurves[count,:]), size(R)
		allCurves[count,:]=R
		count = count+1
	end
	xlabel("Time (minutes)", fontsize = 30)
	ylabel("Clot Amplitude (mm)", fontsize = 30)
	@show size(mean(allCurves, dims=1))
	@show size(T)
	plot(T, vec(mean(allCurves, dims=1)), color="lightblue")
	ax = gca()
	ax[:tick_params]("both",labelsize=24) 
	numSamples = size(solved,1)
	upperlim = mean(allCurves,dims=1)+1.95*std(allCurves,dims=1)/sqrt(numSamples)
	lowerlim = mean(allCurves,dims=1)-1.95*std(allCurves,dims=1)/sqrt(numSamples)
	fill_between(T, vec(upperlim), vec(lowerlim), color = "lightblue", alpha = .5)
	savefig("../figures/solvingSameProbSameCurves_13_03_19.pdf")
end

function analyzeDiffICs(startIdx, endIdx)
	close("all")
	trueICStr = "../solveInverseProb/solveDiffProb/ics_to_match_18_03_19_iter"
	foundICStr ="../solveInverseProb/solveDiffProb/foundIcs_18_03_19_"
	global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	#nominal_ICs=[1106.0, 0.0, 60.0, 0.0, 2686.0, 0.79, 3.95, 0.0, 0.0, 7900.0, 971.7, 1.58, 0.0, 0.0, 0.0, 932.2, 0.4424, 0.0]
	#nominal_experimental=[1.975,15.8,0.553,71.1,134.3,73.47,70.0]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]

	@show size(nominal_experimental), size(nominal_ICs), size(platelets)

	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	errCutoff = 1.0

	exp_condition_names = ["TFPI", "FV", "FVIII", "FIX", "FX", "FXIII", "TAFI", "FXIII"]
	IC_names =["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER","Fraction of Platelets Activated", "FV+FX", "FV+FXa","Prothrombinase Complex", "Fibrin","Plasmin","Fibrinogen", "Plasminogen", "tPA", "uPA", "Fibrin \n monomer", "Protofibril", "antiplasmin", "PAI_1", "Fiber"]
	selidxs = [1,2,3,4,5,6,7,8,9,11,13,14,16,17,20,21,22,23]
	allnames = vcat(exp_condition_names, IC_names)
	numparams = maximum(size(allnames))
#	data = Array{Any,numparams}
#	
#	for j = 1:numparams
#		data[j] = []
#	end
	data = zeros(size(selidxs,1), maximum(size(startIdx:endIdx)))

	for j = startIdx:endIdx
		currFound = readdlm(string(foundICStr, j, ".txt"))
		currGiven = readdlm(string(trueICStr, j, ".txt"),',')
		@show currFound
		@show currGiven
		count = 1
		for k in selidxs
			#@show currFound[k]
			#@show currGiven[k]
			data[count,j]=((currFound[count]/all_nominal[selidxs[count]]-currGiven[count]/all_nominal[selidxs[count]])^2)^.5
			count = count+1
		end
	end
	@show data, size(data)
	meandata =mean(data,dims=1)
#	meandata = Array{Float64,maximum(size(selidxs)), size(data[1],1)}
#	for k =1:maximum(size(selidxs))
#		meandata[k,:]=(data[selidxs[k]])
#	end

	figure(figsize = (20,15))
	#@show size(data[1])
	numSamples = size(data[1],1)
	@show data
	#boxplot(data[selidxs]*100, "k")
	positions = collect(1:maximum(size(selidxs)))
	@show typeof(positions)
	@show typeof(mean(meandata,2)*100)
	@show numSamples
	#errorbars will be std error of the mean
	bar(positions, squeeze(mean(meandata,dims=2),2)*100, yerr = squeeze(std(meandata,dims=2),2)*100/sqrt(numSamples),color = "lightblue")
	ax = gca()
	#ax[:set_yscale]("log")
	ax[:xaxis][:set_ticks](positions)
	#axis([0,numparams,0,3.6])
	ax[:tick_params](labelsize=20)
	#lines and label for kinetic parameters
	ax[:xaxis][:set_ticklabels](allnames[selidxs], rotation = 45, fontsize = 24, horizontalalignment = "right")
	axis("tight")
	#@show size(positions), size(vec(trueICs))
	ax = gca()
	ax[:set_ylim]([0,50])
	ylabel("Initial Concentration\n (% difference of nominal)", fontsize = 36)
	savefig("figs/DiffInICSolvingDiffProb_05_09_18.pdf",bbox_inches="tight")
end
