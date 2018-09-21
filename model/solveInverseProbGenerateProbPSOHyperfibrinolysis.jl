#to solve the inverse problem-we know ROTEM metrics, but can we back out ICs?
using Optim

include("runModel.jl")

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

function objective_six_metrics_weighted(params)
	tic()
	#check for positivity-this matters in unbounded optimization
	if(true in (params .<0))
		calc_obj = 10^8
		return calc_obj
	end
	curr_exp = params[1:8]
	curr_ICs = params[9:end-1]
	curr_platelets = params[end]
	tPA = curr_ICs[12]
	T,R =runModelWithParamsChangeICReturnA(kin_params,curr_ICs,curr_exp,curr_platelets)
	#print(T)
	#plot(T,R)
	#CT,CFT,alpha,MCF,A10,A20,LI30,LI60
	metrics = calculateCommonMetrics(R,T)
	CT = metrics[1]
	CFT = metrics[2]
	alpha = metrics[3]
	MCF = metrics[4]
	#experiment run for 120 mins =2 hours
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
	#@show sel_metrics, calc_obj
	toc()
	return calc_obj

end

function testIfPhysical(R,T)

	#describe the curve we've created'
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	target_AUC = calculateAUC(R,T)

	temp_target=[target_CT, target_CFT, target_alpha, target_MCF,target_MaximumLysis, target_AUC]
	if(true in (temp_target .<0))
		isPhysical = false
	else
		isPhysical = true
	end
	return isPhysical
end

function testIfNominal(R,T)
	#nominal values from Multi-centre investigation on reference ranges for ROTEM thromboelastometry.
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	target_AUC = calculateAUC(R,T)

	nominal = true
	#if(!(target_CT>30/60 && target_CT<246/60))
	#	nominal =false
	#end
	#if(!(target_CFT>40/60 && target_CFT<148/60))
	#	nominal =false
	#end
	if(!(target_MCF> 40 && target_MCF<72))
		nominal = false
	end
	@show target_CT, target_CFT, target_MCF
	return nominal
	


end

function testIfHyperfibrinolysis(R,T)
	#nominal values from Multi-centre investigation on reference ranges for ROTEM thromboelastometry.
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	LI_30 =calculateLysisAtTime(R,T,30.0)
	target_AUC = calculateAUC(R,T)
	LI_90=calculateLysisAtTime(R,T,90.0)
	#as defined by Hyperfibrinolysis, physiologic fibrinolysis, and fibrinolysis shutdown: The spectrum of postinjury fibrinolysis and relevance to antifibrinolytic therapy

	hyperfibrinolysis = false
	if(LI_90<.9 && LI_90>=0.0)
		hyperfibrinolysis = true
	end
	@show LI_90, hyperfibrinolysis
	return hyperfibrinolysis
	


end

function generatePeturbed(lower, upper, cond)
	return cond*lower+rand(size(cond)).*(upper*cond-cond*lower)
end

function runPSOHyperfibrinolysis(seed,iter, numFevals)
	# Load data -
	#global kin_params = readdlm("../parameterEstimation/startingPoint_02_05_18.txt")
	global kin_params = mean(readdlm("../parameterEstimation/best8_02_05_18.txt"),1)
	d = buildCompleteDictFromOneVector(kin_params)
	nominal_ICs = d["INITIAL_CONDITION_VECTOR"]
	#tPA concentration based on https://diapharma.com/tissue-plasminogen-activator-tpa/#Concentration/activity
	nominal_ICs[16]=.0735*10 #2.0#give us some tPA
	d["INITIAL_CONDITION_VECTOR"] = nominal_ICs
	#nominal_ICs=[1106.0, 0.0, 60.0, 0.0, 2686.0, 0.79, 3.95, 0.0, 0.0, 7900.0, 971.7, 1.58, 0.0, 0.0, 0.0, 932.2, 0.4424, 0.0]
	#nominal_experimental=[1.975,15.8,0.553,71.1,134.3,73.47,70.0]
	nominal_experimental = d["FACTOR_LEVEL_VECTOR"]
	platelets = d["PLATELET_PARAMS"][5]

	#@show size(nominal_experimental), size(nominal_ICs), size(platelets)

	#concat together
	all_nominal =vcat(nominal_experimental,nominal_ICs, platelets)
	#create upper and lower bounds
	lbs = zeros(size(all_nominal))
	ups  =2*all_nominal
	#make it possible for things that are zero to move some
	ups[ups.==0]=2.0
	#ups[9]=.1 #limit the amount of FIIa we can start with
	#@show all_nominal[9]
	#bound that we must have some prothrombin
	#lbs[8]=100.0

	#replace zeros in nominal conditions with eps, otherwise NLopt gets upset
	all_nominal[all_nominal.==0]=eps()

	#@show lbs.<ups

	#normal ranges for NATEM from https://www.sciencedirect.com/science/article/pii/S0049384810004081?via%3Dihub
	#CT-461-917 s
	#CFT 105-409 S
	#MCF 26-60 mm
	#alpha 29-65 degrees
	#maximumLysis 8-24 percent
	#maximumVelocity 2-10 mm/min

	#set seed so this is repeatable
	srand(seed)
	#create our data
	#conditions taken from SampleSingleParameterSet adjustments to IC and experimental conditions
	temp_IC = d["INITIAL_CONDITION_VECTOR"]
	temp_exp = d["FACTOR_LEVEL_VECTOR"]
	temp_platelets = d["PLATELET_PARAMS"][5] 
	#with these two stretch factors, we're going between 50% and 150% of normal
	stretchfactorupper = 2.0 #for setting upper and lower bounds of our generated ROTEM curve
	stretchfactorlower = .5	
	genIC =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_IC)
	genExp =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_exp)
	genPlatelets =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_platelets)


	#for hyperfibrinolysis, increase plasmin, plasminogen, tPA, decrease antiplasmin, PAI-1		
	poss_idxs_pos = [13,14,16]
	poss_idxs_neg = [20,21]
	all_poss_idxs =vcat(poss_idxs_pos, poss_idxs_neg)
	#let's only adjust one of these for now
	num_to_adjust =2
	adj_factor =4.0
	for j =1:num_to_adjust
		curr_idx = all_poss_idxs[rand(1:end)]
		if(curr_idx in poss_idxs_pos)
			genIC[curr_idx] = genIC[curr_idx]*adj_factor
		elseif(curr_idx in poss_idxs_neg)
			genIC[curr_idx] = genIC[curr_idx]*1/adj_factor
		end
	end

	genIC = vec(genIC)
	genExp = vec(genExp)
	tPA = genIC[12]
	#@show genIC-temp_IC
	#@show genExp-temp_exp

	T,R =runModelWithParamsChangeICReturnA(kin_params,genIC,genExp,genPlatelets)
	#@show T
	isPhysical = testIfPhysical(R,T)
	isHyperfibrinolysis = testIfHyperfibrinolysis(R,T)
	while(!isPhysical || !isHyperfibrinolysis)
		#this parameters produce a nonsensical results, regenerate and repeat as necceassary 
		println("Checking for feasibility")
		genIC =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_IC)
		genExp =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_exp)
		genPlatelets =generatePeturbed(stretchfactorlower, stretchfactorupper, temp_platelets)
		genIC = vec(genIC)
		genExp = vec(genExp)
		#let time delay range between 400 and 1500
		genMX = 400+rand()*(1500-400)
		tPA = genIC[12]
		T,R =runModelWithParamsChangeICReturnA(kin_params,genIC,genExp,genPlatelets)
		isPhysical = testIfPhysical(R,T)
		isHyperfibrinolysis = testIfHyperfibrinolysis(R,T)
		@show isPhysical, isHyperfibrinolysis, genIC
	end

	#Store our ICS to disk
	println("Suitable conditions found. Writing to disk")
	writedlm(string("../solveInverseProb/ics_to_match_29_07_18_hyperfibrinolysis_iter",iter ,".txt"), hcat(genExp', genIC', genPlatelets), ',')

	#describe the curve we've created'
	metrics = calculateCommonMetrics(R,T)
	target_CT = metrics[1]
	target_CFT = metrics[2]
	target_alpha = metrics[3]
	target_MCF = metrics[4]
	#experiment run for 120 mins =2 hours
	target_MaximumLysis = calculateLysisAtTime(R,T,60.0)
	target_AUC = calculateAUC(R,T)

	#use the midpoint as characteritic scaling
	scale_CT = 500.0
	scale_CFT = 250.0
	scale_alpha = 40.0
	scale_MCF = 40.0
	scale_MaxLysis = .84
	scale_AUC = 2500.0

	global sel_scales = [scale_CT, scale_CFT,scale_alpha, scale_MCF, scale_MaxLysis, scale_AUC]

	#set target-change here

	#if there are more metrics to target, they can be added here 

	global sel_target = [target_CT, target_CFT, target_alpha, target_MCF,target_MaximumLysis, target_AUC]
	@show sel_target
	global weights = [1,1,1.0,1,1,1] 
	#run optimization
#	print("Starting PSO")
#	tic()
	#res = optimize(objective_six_metrics_weighted,lbs,ups, SimulatedAnnealing(), Optim.Options(iterations=numFevals))
#	res = optimize(objective_six_metrics_weighted,lbs,ups, ParticleSwarm(n_particles=40), Optim.Options(iterations=numFevals))
#	toc()
#	print(res)
#	writedlm(string("../solveInverseProb/foundIcs_29_07_18_hyperfibrinolysis_", iter, ".txt"), res.minimizer)
#	writedlm(string("../solveInverseProb/Output_29_07_18_hyperfibrinolysis_", iter, ".txt"),string(res))
end

for p in collect(1:10)
	print("On iter", p, " of 10")
	runPSOHyperfibrinolysis(22+p,p,500)
end

#iterNum =  parse(Int64, ARGS[1])
#runSA(22+iterNum,iterNum,5)




