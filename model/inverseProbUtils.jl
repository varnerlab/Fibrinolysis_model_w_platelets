#utils for analyzing solution from solving inverse problem
using PyPlot
using DelimitedFiles
using Statistics
using Random
using Distributions
using CSV

include("runModel.jl")

function convertFibrinogenTonM(fibrinogen_in_mg_dl::Real)
	mw_fibrinogen = 340.0 #kDa
	kilodaltons_in_a_mg =  6.0221366516752E17
	deci_liters_per_liter =10.0
	NA = 6.0221409e23 #avagadro's constant'
	conv_fibrinogen = fibrinogen_in_mg_dl*1/mw_fibrinogen*kilodaltons_in_a_mg*deci_liters_per_liter*1/NA*10^9
	return conv_fibrinogen
end

function convertPAIToM(PAI1_in_ng_ml::Real)
	mw_PAI1 = 43.0 #kDA
	kilodaltons_in_a_mg =  6.0221366516752E17
	milli_liters_per_liter =1000.0
	NA = 6.0221409e23 #avagadro's constant'
	conv_PAI1 =PAI1_in_ng_ml*1E-6*1/mw_PAI1*kilodaltons_in_a_mg*milli_liters_per_liter*1/NA*10^9
	return conv_PAI1
end

function convertTAFIToM(TAFI_in_ng_ml::Real)
	mw_TAFI = 60 #kDA
	kilodaltons_in_a_mg =  6.0221366516752E17
	milli_liters_per_liter =1000.0
	NA = 6.0221409e23 #avagadro's constant'
	conv_TAFI =TAFI_in_ng_ml*1E-6*1/mw_TAFI*kilodaltons_in_a_mg*milli_liters_per_liter*1/NA*10^9
	return conv_TAFI
end

function convertTFPIToM(TFPI_in_ng_ml::Real)
	mw_TFPI = 32 #kDA
	kilodaltons_in_a_mg =  6.0221366516752E17
	milli_liters_per_liter =1000.0
	NA = 6.0221409e23 #avagadro's constant'
	conv_TFPI =TFPI_in_ng_ml*1E-6*1/mw_TFPI*kilodaltons_in_a_mg*milli_liters_per_liter*1/NA*10^9
	return conv_TFPI
end

function sampleTimeDelay_UW_R()
	mu = 1.66 #minutes
	sigma = sqrt(sqrt(5.87))
	#use truncated normal to prevent negative time delays
	d =Truncated(Normal(mu, sigma), 0, 1000)
	delay = rand(d, 1)[1]
	return delay 
end

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

function testIfPhysicalReturnTrajectory(params,genIC,genExp,genPlatelets)
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

	return isPhysical, temp_target, R
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

function objective_MSE(params)
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	T,R =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	sum = 0
	for j in 1:length(R)
		sum = sum+(R_to_match[j]-R[j])^2
	end
	return sum/length(R)
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

function objective_MSE_moreComplete(params)
	#@show params
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	TFPI = curr_ICs[13]
	TAFI = curr_ICs[14]
	dict = buildCompleteDictFromOneVector(kin_params)
	dict["FACTOR_LEVEL_VECTOR"][1]=TFPI
	dict["FACTOR_LEVEL_VECTOR"][7]=TAFI
	TSTART = 0.0
	TSTOP = 60.0
	step = .1
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	#update our time delay
	est_delay =sampleTimeDelay_UW_R()
	dict["TIME_DELAY"][1]=est_delay
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
	t =sol.t
	X = sol
	#A = convertToROTEM(t,X,tPA)
	tPA = initial_condition_vector[16]
	A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
	sum = 0
#	@show length(A)
#	@show length(R_to_match)
	#plot(t, A)
	for j in 1:length(R_to_match)
		sum = sum+(R_to_match[j]-A[j])^2
	end
	return sum/length(R_to_match)
end

function objective_UW_moreComplete(params::Array)
	#print("here!")
	#@show params
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	TFPI = curr_ICs[13]
	TAFI = curr_ICs[14]
	dict = buildCompleteDictFromOneVector(kin_params)
	dict["FACTOR_LEVEL_VECTOR"][1]=TFPI
	dict["FACTOR_LEVEL_VECTOR"][7]=TAFI
	TSTART = 0.0
	TSTOP = 60.0
	step = .1
	initial_condition_vector = dict["INITIAL_CONDITION_VECTOR"]
	#update our time delay
	est_delay =sampleTimeDelay_UW_R()
	dict["TIME_DELAY"][1]=est_delay
	fbalances(y,p,t)= Balances(t,y,dict) 
	prob = ODEProblem(fbalances, initial_condition_vector, (TSTART,TSTOP))
	sol = solve(prob, alg_hints=[:stiff] , dt = 2.0, dtmax = 1.0, abstol = 1E-6, reltol = 1E-4, force_dtmin=true, saveat = step,maxiters = 1e6)
	t =sol.t
	X = sol
	#A = convertToROTEM(t,X,tPA)
	tPA = initial_condition_vector[16]
	A = convertToROTEMPlateletContribution_UWRescaled(t,X,tPA,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(A,t)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = calculateLysisAtTime(A,t,30.0)
	sel_metrics = [CT,CFT,alpha,MCF,LI30]
	#@show sel_metrics
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
		calc_obj = sqrt(sum)
	end
	#@show calc_obj
	return calc_obj
	
end


function objective_UW(params::Vector, grad::Vector)
	#print("here!")
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	T,R= runModelWithParamsChangeICReturnA_UWRescaled(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = calculateLysisAtTime(R,T,30.0)
	sel_metrics = [CT,CFT,alpha,MCF,LI30]
	#@show sel_metrics
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
		calc_obj = sqrt(sum)
	end
	#@show calc_obj
	return calc_obj
	
end

function objective_UW(params::Array)
	#print("here!")
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	dilution_factor =  0.79
	T,R= runModelWithParamsChangeICReturnA_UWRescaled(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = calculateLysisAtTime(R,T,30.0)
	sel_metrics = [CT,CFT,alpha,MCF,LI30]
	#@show sel_metrics
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
		calc_obj = sqrt(sum)
	end
	#@show calc_obj
	return calc_obj
	
end

function objective_UW_iterative(params::Array)
	#print("here!")
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	dilution_factor =  0.79
	T,R= runModelWithParamsChangeICReturnA_UWRescaled(kin_params,curr_ICs,curr_exp,curr_platelets)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	LI30 = calculateLysisAtTime(R,T,30.0)
	sel_metrics = [MCF,CFT,CT,alpha,LI30]
	#@show sel_metrics[1:metric_number]
	#@show sel_metrics
	#if any of the metrics are negative, something went wrong, penalize our parameters
	if(true in (metrics .<0)) #if any of our metrics are negative
		calc_obj = 1E8
	else
		sum = 0.0
		#use MSE to determine how far away we are from hitting our target
		for j in 1:metric_number
			sum = sum+(sel_metrics[j]/sel_scales[j]-sel_target[j]/sel_scales[j])^2*weights[j]
		end
		#@show params
		calc_obj = sqrt(sum)
	end
	#@show curr_ICs
	#@show sel_metrics[1:metric_number],calc_obj
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
	solved=[collect(1:19);collect(201:206)]
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

function analyzeDiffICs(solvedIdxs)
	close("all")
	#trueICStr = "../solveInverseProb/solveDiffProb_200Evals/ics_to_match_19_03_19_iter"
	#foundICStr ="../solveInverseProb/solveDiffProb_200Evals/foundIcs_19_03_19_"
	trueICStr = "../solveInverseProb/solveDiffProb_300Evals/ics_to_match_23_10_19_iter"
	foundICStr = "../solveInverseProb/solveDiffProb_300Evals/foundIcs_23_10_19_"
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	kin_params=mean(allp, dims=1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	#nominal_ICs=[1106.0, 0.0, 60.0, 0.0, 2686.0, 0.79, 3.95, 0.0, 0.0, 7900.0, 971.7, 1.58, 0.0, 0.0, 0.0, 932.2, 0.4424, 0.0]
	#nominal_experimental=[1.975,15.8,0.553,71.1,134.3,73.47,70.0]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073

	@show size(nominal_experimental), size(nominal_ICs), size(platelets)

	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	errCutoff = 1.0

	exp_condition_names = ["TFPI", "FV", "FVIII", "FIX", "FX", "FXIII", "TAFI", "FXIII"]
	IC_names =["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER","Fraction of Platelets Activated", "FV+FX", "FV+FXa","Prothrombinase Complex", "Fibrin","Plasmin","Fibrinogen", "Plasminogen", "tPA", "uPA", "Fibrin \n monomer", "Protofibril", "antiplasmin", "PAI_1", "Fiber"]
	selidxs = [1,3,4,6,7,8,9,11,13,14,17,22,23,24]
	allnames = vcat(exp_condition_names, IC_names)
	@show allnames[selidxs]
	numparams = maximum(size(allnames))
#	data = Array{Any,numparams}
#	
#	for j = 1:numparams
#		data[j] = []
#	end
	data = zeros(size(selidxs,1), maximum(size(solvedIdxs)))
	@show size(data)

	itercount = 1
	for j in solvedIdxs
		currFound = readdlm(string(foundICStr, j, ".txt"))
		currGiven = readdlm(string(trueICStr, j, ".txt"),',')
		@show currFound
		@show currGiven
		count = 1
		for k in selidxs
			#@show currFound[k]
			#@show currGiven[k]
			data[count,itercount]=((currFound[k]/all_nominal[k]-currGiven[k]/all_nominal[k])^2)^.5
			count = count+1
		end
		itercount = itercount+1
	end
	@show data, size(data)
	meandata =mean(data,dims=1)
#	meandata = Array{Float64,maximum(size(selidxs)), size(data[1],1)}
#	for k =1:maximum(size(selidxs))
#		meandata[k,:]=(data[selidxs[k]])
#	end

	figure(figsize = (20,15))
	#@show size(data[1])
	numSamples = size(data,2)
	@show data
	#boxplot(data[selidxs]*100, "k")
	positions = collect(1:maximum(size(selidxs)))
	#errorbars will be std error of the mean
	std_err=vec(std(data,dims=2)*100/sqrt(numSamples))
	@show size(positions)
	@show size(std_err)
	@show size(vec(mean(data,dims=2)*100))
	bar(positions, vec(mean(data,dims=2)*100), yerr = std_err,color = "lightblue")
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
	ax[:set_ylim]([0,150])
	ylabel("Initial Concentration\n (% difference of nominal)", fontsize = 36)
	savefig("../figures/DiffInICSolvingDiffProb_23_10_19_N52.pdf",bbox_inches="tight")
end

function analyzeDiffICs_UW()
	close("all")
	date_str = "Acc_11_11_19"
	allp = readdlm("../LOOCV/bestparamsForBatch_10_14_02_19.txt")
	#global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	kin_params=mean(allp, dims=1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	#nominal_ICs=[1106.0, 0.0, 60.0, 0.0, 2686.0, 0.79, 3.95, 0.0, 0.0, 7900.0, 971.7, 1.58, 0.0, 0.0, 0.0, 932.2, 0.4424, 0.0]
	#nominal_experimental=[1.975,15.8,0.553,71.1,134.3,73.47,70.0]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]
	#add some tPA
	nominal_ICs[16]=.073

	@show size(nominal_experimental), size(nominal_ICs), size(platelets)

	#concat together
	all_nominal =vcat(nominal_ICs, nominal_experimental,platelets)
	errCutoff = 1.0

	exp_condition_names = ["TFPI", "FV", "FVIII", "FIX", "FX", "FXIII", "TAFI", "FXIII"]
	IC_names =["FII", "FIIa", "PC", "APC", "ATIII", "TM", "TRIGGER","Fraction of Platelets Activated", "FV_FX", "FV_FXa","Prothrombinase Complex", "Fibrin","Plasmin","Fibrinogen", "plasminogen", "tPA", "uPA", "Fibrin \n monomer", "Protofibril", "antiplasmin", "PAI-1", "Fiber"]
	allnames = vcat(IC_names,exp_condition_names)
	numparams = maximum(size(allnames))

	#get number of sims by counting files
	all_fns=readdir("../solveInverseProb/UW/")
	num_sims =sum(occursin.(date_str,all_fns))
	data = zeros(length(allnames), num_sims)
	@show size(data)
	#which indexes correspond to the initial conditions we've been solving for
	condition_idxs =[1,3,5,14,15,9,21,23]

	pathToData = "/home/rachel/Documents/washington_burns/PatientsWTEGAndMoreFactors.csv"
	all_data = CSV.read(pathToData)
	for j in 1:size(all_data, 1)
		curr_id = all_data.id[j]
		curr_hr = all_data.hour[j]
		curr_res =  readtable(string("../solveInverseProb/UW/",date_str,"_comparisonPSOPatient",curr_id, "Hour", curr_hr,".txt"), separator='\t')
		for k in 1:size(curr_res)[1]
			condition_idx =findfirst(x->x==string(curr_res[!,:Name][k]), allnames)
			#@show curr_res[k,:Actual],curr_res[k,:Estimated], condition_idx
			data[condition_idx,j]=((curr_res[k,:Actual]/all_nominal[condition_idx]-curr_res[k,:Estimated]/all_nominal[condition_idx])^2)^.5
		end
	end
	@show data, size(data)
	meandata =mean(data,dims=1)

	figure(figsize = (20,15))
	#@show size(data[1])
	numSamples = size(data,2)
	@show data[condition_idxs, :]
	#boxplot(data[selidxs]*100, "k")
	positions = collect(1:maximum(size(condition_idxs)))
	#errorbars will be std error of the mean
	std_err=vec(std(data,dims=2)*100/sqrt(numSamples))
	@show size(positions)
	@show size(std_err)
	@show size(vec(mean(data,dims=2)*100))
	bar(positions, vec(mean(data[condition_idxs],dims=2)*100), yerr = std_err[condition_idxs],color = "lightblue")
	ax = gca()
	#ax[:set_yscale]("log")
	ax[:xaxis][:set_ticks](positions)
	#axis([0,numparams,0,3.6])
	ax[:tick_params](labelsize=20)
	#lines and label for kinetic parameters
	ax[:xaxis][:set_ticklabels](allnames[condition_idxs], rotation = 45, fontsize = 24, horizontalalignment = "right")
	axis("tight")
	#@show size(positions), size(vec(trueICs))
	ax = gca()
	ax[:set_ylim]([0,300])
	ylabel("Initial Concentration\n (% difference of nominal)", fontsize = 36)
	savefig(string("../figures/DiffInICSolvingDiffProb",date_str,".pdf"),bbox_inches="tight")
end
